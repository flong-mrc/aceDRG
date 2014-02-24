
/* 
 * File:   MolSdfFile.cpp
 * Author: flong
 *
 * Created on November 01, 2012, 11:25 AM
 */

#include "MolSdfFile.h"

namespace LIBMOL
{
    MolSdfFile::MolSdfFile():  
            molName(NullString),
            creatProg(NullString),
            hasCoords(true),
            hasConnect(false),
            hasH(false),
            containMetal(false),
            itsCurMol(NullPoint)
    {
    }
    
    MolSdfFile::MolSdfFile(Name tFname, std::ios_base::openmode tOpenMode) :
              molName(NullString),
              creatProg(NullString),
              hasCoords(true),
              hasConnect(false),
              hasH(false),
              containMetal(false),
              itsCurMol(NullPoint)
    {
        if (tOpenMode == std::ios::in)
        {
            inFile.open(tFname.c_str(), tOpenMode);
            setupSystem();
        }
        else
        {
            outFile.open(tFname.c_str(), tOpenMode);
        }    
    }
    
    MolSdfFile::MolSdfFile(FileName tFname, std::ios_base::openmode tOpenMode) :
              molName(NullString),
              creatProg(NullString),
              hasCoords(false),
              hasConnect(false),
              hasH(false),
              itsCurMol(NullPoint)
    {
        
        
        if (tOpenMode == std::ios::in)
        {
            
            inFile.open(tFname, tOpenMode);
            setupSystem();
        }
        else
        {
            outFile.open(tFname, tOpenMode);
        }    
    }
    
    MolSdfFile::~MolSdfFile()
    {
        deleteCurMol();
    }
    
    void MolSdfFile::setupSystem()
    {
        
        if (inFile.is_open() )
        {  
            // related to records
            bool  lStart=true, lEnd=false;
            ID    curProp;
            int   nAtoms=0, nBonds=0;
            int   iLineIdx =0, tAtomIdx =0, tBondIdx=0;
            std::map<ID, std::vector<int> > tIDs;
            
            while(!inFile.eof())
            {       
                std::string tRecord="";
                std::getline(inFile, tRecord);
                tRecord = TrimSpaces(tRecord);
                std::cout << tRecord << std::endl;
                
                if (lStart)
                {
                    createCurMol(); 
                    iLineIdx =0;
                    tAtomIdx = 0;
                    tBondIdx = 0;
                    lStart   = false;
                    
                }
                
                std::size_t tFind1 = tRecord.find("$$$$");
                std::size_t tFind2 = tRecord.find("END");
                std::size_t tFind3 = tRecord.find(">");
                
                if (tFind1==std::string::npos)
                {
                    if(lEnd)
                    {
                        if (tFind3 ==std::string::npos && tRecord.size() !=0)
                        {
                            itsCurMol->propData[curProp].push_back(tRecord);
                        }
                        else
                        { 
                            std::vector<std::string> tKeyBuf;
                            StrTokenize(tRecord, tKeyBuf, '<');
                            if (tKeyBuf.size() ==2)
                            {
                                curProp = tKeyBuf[1].erase(tKeyBuf[1].size()-1, 1);
                            }
                            else
                            {
                                std::cout << "The property key word in SDF file is "
                                          << tRecord << std::endl;
                            }
                        }
                    }
                    else if(tFind2 !=std::string::npos)
                    {
                        iLineIdx =0;
                        lEnd = true;
                    }
                    else
                    {
                        if (iLineIdx <=2)
                        {
                            itsCurMol->comments.push_back(tRecord);
                        }
                        else if (iLineIdx ==3)
                        {
                            std::vector<std::string> tBuf;
                            StrTokenize(tRecord, tBuf);
                            nAtoms = StrToInt(tBuf[0]);
                            nBonds = StrToInt(tBuf[1]);
                        }
                        else if (iLineIdx > 3)
                        {
                            std::vector<std::string> tBuf;
                            StrTokenize(tRecord, tBuf);
                            if (iLineIdx <= nAtoms+3)
                            {
                                AtomDict aAtom;
                                aAtom.seriNum = tAtomIdx;
                                std::cout << "aAtom.seriNum " << tAtomIdx << std::endl;
                                aAtom.coords[0] = StrToReal(tBuf[0]);
                                aAtom.coords[1] = StrToReal(tBuf[1]);
                                aAtom.coords[2] = StrToReal(tBuf[2]);
                                aAtom.chemType = TrimSpaces(tBuf[3]);
                                tIDs[aAtom.chemType].push_back(aAtom.seriNum);
                                int tSize = (int)tIDs[aAtom.chemType].size();
                                aAtom.id       = aAtom.chemType + IntToStr(tSize);
                                itsCurMol->atoms.push_back(aAtom);
                                tAtomIdx++;
                            }
                            else if (iLineIdx > (nAtoms+3) && iLineIdx <=(nAtoms+nBonds+3))
                            {
                                BondDict aBond;
                                int t1 = StrToInt(tBuf[0])-1;
                                int t2 = StrToInt(tBuf[1])-1;
                                aBond.atomsIdx.push_back(t1);
                                aBond.atomsIdx.push_back(t2);
                                aBond.atoms.push_back(itsCurMol->atoms[t1].id);
                                aBond.atoms.push_back(itsCurMol->atoms[t2].id);
                                itsCurMol->atoms[t1].connAtoms.push_back(t2);
                                itsCurMol->atoms[t2].connAtoms.push_back(t1);
                                aBond.fullAtoms[itsCurMol->atoms[t1].id] = t1;
                                aBond.fullAtoms[itsCurMol->atoms[t2].id] = t2;
                                aBond.order   = tBuf[2];
                                aBond.seriNum = tBondIdx;
                                std::cout << "aBond.seriNum " << aBond.seriNum << std::endl;
                                itsCurMol->bonds.push_back(aBond);
                                tBondIdx++;
                            }
                        }
                        
                        iLineIdx++;
                    }
                    
                }
                else
                {
                    allMols.push_back(*itsCurMol);
                    deleteCurMol();
                    nAtoms=0; 
                    nBonds=0;
                    iLineIdx =0;
                    tAtomIdx =0; 
                    tBondIdx =0;
                    lStart   =true;
                    lEnd     =false;
                }
            }
            if (itsCurMol->atoms.size() !=0)
            {
                allMols.push_back(*itsCurMol);
            }
            inFile.close();
        }
       
        for (unsigned i=0; i < allMols.size(); i++)
        {
            addHAtomToMols(i);
            setAtomsBondingAndChiralCenter(i);
            setChiral(i);
        }
       
        
        // check 
        std::cout << "There are " << (int)allMols.size() 
                << " molecule(s) in the SDF file " << std::endl;
        int i = 0;
        for (std::vector<Molecule>::iterator iM=allMols.begin();
                iM !=allMols.end(); iM++)
        {
            i++;
            std::cout << "For Molecule No. " << i << " : " << std::endl;
            std::cout << "Its header information is " << std::endl;
            for (std::vector<std::string>::iterator iComm=iM->comments.begin();
                    iComm != iM->comments.end(); iComm++)
            {
                std::cout<< *iComm << std::endl;
            }
            
            std::cout << "The molecule has " <<  iM->atoms.size() 
                      << " atoms, including " << iM->extraHAtoms.size()
                      << " newly added H atoms. "<< std::endl;
            std::cout << "These atoms are : " << std::endl;
           
            for (std::vector<AtomDict>::iterator iA=iM->atoms.begin();
                    iA !=iM->atoms.end(); iA++)
            {
                std::cout << "atom " << iA->seriNum +1  << "\t" << iA->chemType 
                          << " id : " << iA->id << "\t"
                        << "\t" << iA->coords[0] << "\t" 
                        << iA->coords[1] << "\t"
                        << iA->coords[2] << std::endl;
                if (iA->connAtoms.size() !=0)
                {
                    std::cout << "It connects to the following atoms: " <<std::endl;
                    for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                             iNB !=iA->connAtoms.end(); iNB++)
                    {
                        std::cout << "Atom " << *iNB+1 << std::endl;
                    }
                }
            }
            std::cout << "There are " << (int)iM->bonds.size() 
                    << " bonds in the molecule and they are :" << std::endl;
            for (std::vector<BondDict>::iterator iB=iM->bonds.begin();
                    iB != iM->bonds.end(); iB++)
            {
                std::cout << "bond " << iB->seriNum + 1 << std::endl
                          << " between atom " << iB->atomsIdx[0]+1
                          << " and atom " << iB->atomsIdx[1]+1
                          << " of order " << iB->order << std::endl;
            }
            
            
            if (iM->propData.size() !=0)
            {
                std::cout << "The molecule has following properties: " << std::endl;
                for (std::map<ID, std::vector<ID> >::iterator iProp=iM->propData.begin();
                        iProp !=iM->propData.end(); iProp++)
                {
                    std::cout << "<" << iProp->first << ">" << std::endl;
                    for (std::vector<ID>::iterator iVal=iProp->second.begin();
                            iVal !=iProp->second.end(); iVal++)
                    {
                        std::cout << *iVal << std::endl;
                    }
                }
            }
        }
        
    }
    
    void MolSdfFile::createCurMol()
    {
        deleteCurMol();
        itsCurMol = new Molecule();
    }
            
    void MolSdfFile::deleteCurMol()
    {
        if (itsCurMol != NullPoint)
        {
            delete itsCurMol;
            itsCurMol = NullPoint;
        }
    }
    
    void MolSdfFile::addHAtomToMols(int  tIdxMol)
    {
        
        PeriodicTable aPTab;
        std::vector<std::string> orgElems;
        orgElems.push_back("C");
        orgElems.push_back("O");
        orgElems.push_back("N");
        orgElems.push_back("S");
        
        for (std::vector<AtomDict>::iterator iA=allMols[tIdxMol].atoms.begin();
                iA !=allMols[tIdxMol].atoms.end(); iA++)
        {
            
            if (std::find(orgElems.begin(), orgElems.end(), iA->chemType)
                 !=orgElems.end())
            {
                // need to check protonated form
                REAL addH =checkProtonated(iA, tIdxMol);
                if (addH >0)
                {
                    addHAtoms(tIdxMol, iA->seriNum, addH);
                }
            }
            
            /*
            if (aPTab.elements.find(iA->chemType) !=aPTab.elements.end())
            {
                if (tVal !=aPTab.elements[iA->chemType]["val"] 
                     && std::find(aPTab.extraValences[iA->chemType].begin(),
                                  aPTab.extraValences[iA->chemType].end(), tVal)
                                  ==aPTab.extraValences[iA->chemType].end())
                {
                    if (tVal < aPTab.elements[iA->chemType]["val"])
                    {
                        REAL diffV= aPTab.elements[iA->chemType]["val"]-tVal;
                        addHAtoms(tIdxMol, iA->seriNum, diffV);
                    }
                    else if (tVal > aPTab.elements[iA->chemType]["val"])
                    {
                        std::cout << "Atom " << iA->seriNum+1 << ", " 
                                  << iA->chemType << std::endl;
                        std::cout << "Its valence by bonding order " << tVal 
                                  << ", its valence in Periodic table " 
                                  << aPTab.elements[iA->chemType]["val"] << std::endl;
                        std::cout << "Warning: the number of bonds binding atom " 
                                  << iA->seriNum+1 << " is larger than its valence "
                                  << std::endl;
                    }
                }
            }
            else
            {
                std::cout << "Element type " << iA->chemType 
                          << " can not be found in the Periodic Table"
                          << std::endl;
                exit(1);
            }
             */
        }
        
        // Now add newly generated H atoms to atoms in the molecule
        // and add associated bonds to the bond set in the molecule
        for (std::vector<AtomDict>::iterator iHA=allMols[tIdxMol].extraHAtoms.begin();
                iHA !=allMols[tIdxMol].extraHAtoms.end(); iHA++)
        {
            allMols[tIdxMol].atoms.push_back(*iHA);
            BondDict aBo;
            aBo.atomsIdx.push_back(iHA->connAtoms[0]);
            aBo.atomsIdx.push_back(iHA->seriNum);
            aBo.order = "1";
            aBo.seriNum = (int)allMols[tIdxMol].bonds.size();
            aBo.atoms.push_back(allMols[tIdxMol].atoms[iHA->connAtoms[0]].id);
            aBo.atoms.push_back(iHA->id);
            aBo.fullAtoms[iHA->id] = iHA->seriNum;
            aBo.fullAtoms[allMols[tIdxMol].atoms[iHA->connAtoms[0]].id] = allMols[tIdxMol].atoms[iHA->connAtoms[0]].seriNum;
            allMols[tIdxMol].bonds.push_back(aBo);
        }
    }
    
    REAL MolSdfFile::checkProtonated(std::vector<AtomDict>::iterator tAt, 
                                     int                  tMolIdx)
    {
        // consider protonated state based atom's property
        REAL addedH =0.0;
        if (tAt->chemType.compare("O")==0)
        {
            addedH=checkProtonateO(tAt, allMols[tMolIdx]);
        }
        else if (tAt->chemType.compare("S")==0)
        {
            addedH=checkProtonateS(tAt, allMols[tMolIdx]);
        }
        else if (tAt->chemType.compare("N")==0)
        {
            addedH=checkProtonateN(tAt, allMols[tMolIdx]);
        }
        else if (tAt->chemType.compare("C")==0)
        {
            addedH=checkProtonateC(tAt, allMols[tMolIdx]);
        }
        
        return addedH;
    }
    
    
    REAL MolSdfFile::checkProtonated(std::vector<AtomDict>::iterator tAt, 
                                     int                  tMolIdx,
                                     REAL                 tTotalVal, 
                                     REAL tPka,           REAL tPh)
    {
        // consider protonated state based atom's property
        REAL addedH=0.0;
        
        if (tAt->chemType.compare("O")==0)
        {
            addedH=checkProtonateO(tAt, allMols[tMolIdx], tPka, tPh);
        }
        else if (tAt->chemType.compare("S")==0)
        {
            addedH=checkProtonateS(tAt, allMols[tMolIdx], tPka, tPh);
        }
        else if (tAt->chemType.compare("N")==0)
        {
            addedH=checkProtonateN(tAt, allMols[tMolIdx], tPka, tPh);
        }
        else if (tAt->chemType.compare("C")==0)
        {
            addedH=checkProtonateN(tAt, allMols[tMolIdx], tPka, tPh);
        }
        
        return addedH;
    }
    
    /*
    REAL MolSdfFile::getBondOrder(int tIdxMol, int tIdx1, int tIdx2)
    {
        REAL tOrd = -1.0;
        
        for (std::vector<BondDict>::iterator iB=allMols[tIdxMol].bonds.begin();
                iB !=allMols[tIdxMol].bonds.end(); iB++)
        {
            if ((iB->atomsIdx[0]==tIdx1 && iB->atomsIdx[1]==tIdx2)
                || (iB->atomsIdx[0]==tIdx2 && iB->atomsIdx[1]==tIdx1))
            {
                tOrd = StrToReal(iB->order);
                if (tOrd ==4.0)
                {
                    tOrd = 1.5;
                }
                
                std::cout << "Bond order " << iB->order << std::endl
                          << tOrd << std::endl;
                break;
            }
        }
        
        return tOrd;
    }
    */
    
    void MolSdfFile::addHAtoms(int tIdxMol, int tIdxAtm, REAL tNumH)
    {
        for (int i=0; i < (int)tNumH; i++)
        {
            AtomDict aH;
            aH.chemType = "H";
            aH.id  = aH.chemType + allMols[tIdxMol].atoms[tIdxAtm].id+IntToStr(i+1);
            aH.seriNum =  (int)allMols[tIdxMol].atoms.size()
                         +(int)allMols[tIdxMol].extraHAtoms.size();
            aH.connAtoms.push_back(allMols[tIdxMol].atoms[tIdxAtm].seriNum);
            allMols[tIdxMol].atoms[tIdxAtm].connAtoms.push_back(aH.seriNum);
            allMols[tIdxMol].atoms[tIdxAtm].connHAtoms.push_back(aH.seriNum);
            
            allMols[tIdxMol].extraHAtoms.push_back(aH);
        }
        allMols[tIdxMol].hasCoords = false;
        
        std::cout << tNumH << " H atoms have been added to bind atom " 
                  << allMols[tIdxMol].atoms[tIdxAtm].id << std::endl;  
    }
    
    int MolSdfFile::getNumOxyConnect(int tIdxMol, std::vector<AtomDict>::iterator iA)
    {
        int nO=0;
        for (std::vector<int>::iterator iC=iA->connAtoms.begin();
                iC !=iA->connAtoms.end(); iC++)
        {
            if (allMols[tIdxMol].atoms[*iC].chemType.compare("O")==0)
            {
                nO++;
            }
        }
        return nO;
    }
    
    void MolSdfFile::setAtomsBondingAndChiralCenter(int tIdxMol)
    {
                // First round
        for (std::vector<AtomDict>::iterator iAt = allMols[tIdxMol].atoms.begin();
                iAt != allMols[tIdxMol].atoms.end(); iAt++)
        {
            int t_len =0;
            for (std::vector<int>::iterator iConn=iAt->connAtoms.begin();
                    iConn !=iAt->connAtoms.end(); iConn++)
            {
                if(!allMols[tIdxMol].atoms[*iConn].isMetal)
                {
                    t_len++;
                }
            }

            if (iAt->chemType.compare("C")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    if (iAt->chiralIdx ==0)
                    {
                        iAt->chiralIdx  = 2;
                    }
                    iAt->bondingIdx = 3;
                }
                else if (t_len ==3)
                {
                    iAt->chiralIdx  = 0;
                    iAt->bondingIdx = 2;
                }
                else if(t_len==2)
                {
                    iAt->chiralIdx  = 0;
                    if (getNumOxyConnect(tIdxMol, iAt)==1)
                    {
                        // water is removed 
                        iAt->bondingIdx=2;
                    }
                    else
                    {
                        iAt->bondingIdx=1;
                    }
                }
            }
            else if (iAt->chemType.compare("N")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    if (iAt->chiralIdx ==0)
                    {
                        iAt->chiralIdx  = 2;
                    } 
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3;  
                }
                //else if (t_len==3) // temp
                //{
                //    iAt->chiralIdx  = -1;
                //    iAt->bondingIdx =  2;
                //}
                else if (t_len ==2)
                {
                    iAt->chiralIdx  = 0;
                    iAt->bondingIdx = 2;
                } 
                else if (t_len==1)
                {
                    // triple bond 
                    iAt->chiralIdx = 0;
                    iAt->bondingIdx= 1;
                }
            }
            else if (iAt->chemType.compare("B")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3;
                }
            }
            else if (iAt->chemType.compare("O")==0)
            {
                if ((int)iAt->connAtoms.size()==2)
                {
                    iAt->bondingIdx = 2;
                }
            }
            else if (iAt->chemType.compare("SI")==0 
                    || iAt->chemType.compare("P")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    
                    if (iAt->chiralIdx ==0)
                    {
                        std::vector<ID> atps;
                        for (std::vector<int>::iterator iNA=iAt->connAtoms.begin();
                                iNA !=iAt->connAtoms.end(); iNA++)
                        {
                            if (std::find(atps.begin(), atps.end(), allMols[tIdxMol].atoms[*iNA].chemType)==atps.end())
                            {
                                atps.push_back(allMols[tIdxMol].atoms[*iNA].chemType);
                            }
                        }
                        if ((int)atps.size() >2)
                        {
                            iAt->chiralIdx  = 2;
                        }
                        else
                        {
                            iAt->chiralIdx =0;
                        }
                        // iAt->chiralIdx  = 2;
                    }
                   
                    iAt->bondingIdx = 3; 
                }
                else if (t_len==3)
                {
                    if (iAt->chiralIdx ==0)
                    {
                        iAt->chiralIdx  = 2;
                    }
                    
                    iAt->bondingIdx = 2; 
                }
            }
            else if (iAt->chemType.compare("S")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4 || t_len==3)
                {
                    if (iAt->chiralIdx ==0)
                    {
                        iAt->chiralIdx  = 2;
                    }
                    iAt->bondingIdx = 3; 
                }
            }
        }
        
        // more conditions 
        
        for (std::vector<AtomDict>::iterator iAt = allMols[tIdxMol].atoms.begin();
                iAt != allMols[tIdxMol].atoms.end(); iAt++)
        {   
            int t_len =0;
            for (std::vector<int>::iterator iConn=iAt->connAtoms.begin();
                    iConn !=iAt->connAtoms.end(); iConn++)
            {
                if(!allMols[tIdxMol].atoms[*iConn].isMetal)
                {
                    t_len++;
                }
            }

            if (iAt->chemType.compare("N")==0 || iAt->chemType.compare("B")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==3)
                {
                    bool l_sp2 = false;
                    for (std::vector<int>::iterator iCA=iAt->connAtoms.begin();
                            iCA != iAt->connAtoms.end(); iCA++)
                    {
                        if(allMols[tIdxMol].atoms[*iCA].bondingIdx == 2)
                        {
                            l_sp2 = true;
                        }
                    }
                    if (l_sp2)
                    {
                        // Now we can say this atom is in sp2 orbits 
                        iAt->chiralIdx  =  0;
                        iAt->bondingIdx =  2;
                    }
                    else
                    {
                        if (iAt->chiralIdx ==0)
                        {
                            iAt->chiralIdx  = 2;
                        }
           
                        iAt->bondingIdx =  3;
                    }
                } 
            }
        }
        
        // Further check if a chiral center is a real one
        for (std::vector<AtomDict>::iterator iA=allMols[tIdxMol].atoms.begin();
                iA != allMols[tIdxMol].atoms.end(); iA++)
        {
            if (iA->chiralIdx !=0)
            {
                
                std::vector<ID> chirRAtms;
                for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                        iNB != iA->connAtoms.end(); iNB++)
                {
                    std::size_t tFind = allMols[tIdxMol].atoms[*iNB].chemType.find("H");
                    if (tFind !=std::string::npos)
                    {
                        chirRAtms.push_back(allMols[tIdxMol].atoms[*iNB].id);
                    }
                }
                
                if ((int)chirRAtms.size() >1 && (int)iA->connAtoms.size() <=4)
                {
                    iA->chiralIdx = 0;
                }
            }
        }
        
        /*
                // First round
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {
            int t_len =0;
            for (std::vector<int>::iterator iConn=iAt->connAtoms.begin();
                    iConn !=iAt->connAtoms.end(); iConn++)
            {
                if(!allAtoms[*iConn].isMetal)
                {
                    t_len++;
                }
            }
            if (iAt->chemType.compare("C")==0)
            {
                //int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3;
                }
                else if (t_len ==3)
                {
                    iAt->chiralIdx  = -1;
                    iAt->bondingIdx = 2;
                } 
            }
            else if (iAt->chemType.compare("N")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3;  
                }
                //else if (t_len==3) // temp 
                //{   // should do on the next round when all NB atoms are set
                //    iAt->chiralIdx  = -1;
                //    iAt->bondingIdx =  2;
                // }
                else if (t_len ==2)
                {
                    iAt->chiralIdx  = -1;
                    iAt->bondingIdx =  2;
                } 
            }
            else if (iAt->chemType.compare("B")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3;
                }
            }
            else if (iAt->chemType.compare("SI")==0 
                    || iAt->chemType.compare("P")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3; 
                }
                else if (t_len==3)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 2; 
                }
            }
            else if (iAt->chemType.compare("S")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4 || t_len==3)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3; 
                }
            }
        }
        
        // more conditions 
        
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {            
            if (iAt->chemType.compare("N")==0 || iAt->chemType.compare("B")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                int t_len =0;
                for (std::vector<int>::iterator iConn=iAt->connAtoms.begin();
                     iConn !=iAt->connAtoms.end(); iConn++)
                {
                    if(!allAtoms[*iConn].isMetal)
                    {
                        t_len++;
                    }
                }
                if(t_len==3)
                {
                    
                    bool l_sp2 = false;
                    for (std::vector<int>::iterator iCA=iAt->connAtoms.begin();
                            iCA != iAt->connAtoms.end(); iCA++)
                    {
                        if(allAtoms[*iCA].bondingIdx == 2)
                        {
                            l_sp2 = true;
                        }
                    }
                    if (l_sp2)
                    {
                        // Now we can say this atom is in sp2 orbits 
                        iAt->chiralIdx  = -1;
                        iAt->bondingIdx =  2;
                    }
                    else
                    {
                        iAt->chiralIdx  =  1;
                        iAt->bondingIdx =  3;
                    }
                } 
            }
        }
 
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA != allAtoms.end(); iA++)
        {
            if (iA->chiralIdx ==1)
            {
                std::vector<ID> chirRAtms;
                for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                        iNB != iA->connAtoms.end(); iNB++)
                {
                    std::size_t tFind = allAtoms[*iNB].chemType.find("H");
                    if (tFind !=std::string::npos)
                    {
                        chirRAtms.push_back(allAtoms[*iNB].id);
                    }
                }
                if ((int)chirRAtms.size() >1 && (int)iA->connAtoms.size() <=4)
                {
                    iA->chiralIdx = 0;
                }
            }
        }
        */
        // No need for the third round, those could be defined in 
        // the first round
        // Check
        /*
        std::cout << "Chiral and plane feather for atoms in the system" 
                  << std::endl;
        
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {
            if (iAt->chiralIdx == -1)
            {
                std::cout << "Atom " << iAt->id << " may be in planes " 
                        << std::endl;
            }
            else if (iAt->chiralIdx == 1)
            {
                std::cout << "Atom " << iAt->id << " may be in a chiral center "
                        << std::endl;
            }
            else if (iAt->chiralIdx==2)
            {
                std::cout << "Atom " << iAt->id 
                        << " may be in a chiral center but the volume sign undefined "
                        << std::endl;
            }
            else 
            {
                std::cout << "Atom " << iAt->id << " is not a chiral center "
                        << std::endl;
            }
        } 
         */  
        
    }
    
    void MolSdfFile::setChiral(int tIdxMol)
    {
        for (std::vector<AtomDict>::iterator iA=allMols[tIdxMol].atoms.begin();
                iA !=allMols[tIdxMol].atoms.end(); iA++)
        {
            if (iA->chiralIdx !=0)
            {
                int tChIdx = (int)allMols[tIdxMol].chirals.size();
                iA->inChirals.push_back(tChIdx);
                ChiralDict aCh;
                aCh.atoms.push_back(iA->seriNum);
                for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                        iNB !=iA->connAtoms.end(); iNB++)
                {
                    aCh.atoms.push_back(*iNB);
                    allMols[tIdxMol].atoms[*iNB].inChirals.push_back(tChIdx);
                }
                aCh.setMutTable(iA->chiralIdx);
                allMols[tIdxMol].chirals.push_back(aCh);
            }
        }
    }
    
    
    extern REAL checkProtonateO(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol)
    {
        REAL aNumH = 0.0;
        StrUpper(tIA->resName);
        // Currently using PH=7 
        // 1. O in Carboxy-Terminus, PKa=2, in ASP and GLU pKa=4
        // No H need to be added
        
        // 2. O in Tyrosine, Pka=10
        if (tIA->connAtoms.size()==1 && tIA->charge==0 )
        {
            if(tIA->resName.compare("TYR")==0)
            {
                aNumH=1.0;
            }
        }
        
        return aNumH;
    }
    
    extern REAL checkProtonateO(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol,
                                REAL tPka,           REAL tPh)
    {
        REAL aNumH=0.0;
        
        return aNumH;
    }
    
    extern REAL checkProtonateN(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol)
    {
        REAL aNumH = 0.0;
        StrUpper(tIA->resName);
        REAL tVal = getTotalBondOrder(tMol, tIA);
        // std::cout << "tVal for N " << tVal << std::endl;
        
        if (tVal < 4)
        {
            if(tIA->resName.compare("LYS")==0 
               || tIA->resName.compare("ARG")==0)
            {
                REAL tD = 4 - tVal;  
                aNumH=tD;
                tIA->formalCharge = tD;
            }
            else
            {
                // check Amino-Terminus, pKa=10 
                int tC=0, tH=0, tOther=0;
                for (std::vector<int>::iterator iNB=tIA->connAtoms.begin();
                        iNB != tIA->connAtoms.end(); iNB++)
                {
                    if(tMol.atoms[*iNB].chemType.compare("C")==0)
                    {
                        tC++;
                    }
                    else if(tMol.atoms[*iNB].chemType.compare("H")==0)
                    {
                        tH++;
                    }
                    else
                    {
                        tOther++;
                    }
                }
                
                if (tC==1 && tOther==0)
                {
                    //It is Amino-Terminus, need to put a charge here
                    aNumH =4 - tVal;
                    
                }
                
            }
        }
        
        return aNumH;
    }
    
    extern REAL checkProtonateN(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol,
                                REAL tPka,           REAL tPh)
    {
        REAL aNumH=0.0;
        
        return aNumH;
    }
    
    extern REAL checkProtonateS(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol)
    {
        REAL aNumH = 0.0;
        
        REAL tVal = getTotalBondOrder(tMol, tIA);
        
        if (tVal==1 && tIA->formalCharge==0 )
        {
            aNumH=1.0;
        }
        
        return aNumH;
    }
    
    extern REAL checkProtonateS(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol,
                                REAL tPka,           REAL tPh)
    {
        REAL aNumH=0.0;
        
        return aNumH;
    }
    
    extern REAL checkProtonateC(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol)
    {
        REAL aNumH = 0.0;
        
        REAL tVal = getTotalBondOrder(tMol, tIA);
        
        if (tVal < 4)
        {
            aNumH = 4-tVal;
        }
        else if (tVal > 4)
        {
            std::cout << "Warning: Atom " << tIA->id << " has bond order "
                      << tVal << std::endl;
        }
        
        return aNumH;
    }
    
    extern REAL checkProtonateC(std::vector<AtomDict>::iterator tIA, 
                                Molecule   & tMol,
                                REAL tPka,           REAL tPh)
    {
        REAL aNumH=0.0;
        
        return aNumH;
        
    }
    
    extern REAL getTotalBondOrder(Molecule   & tMol, 
                                  std::vector<AtomDict>::iterator tIA)
    {
        REAL tVal = 0.0;
        for (std::vector<int>::iterator iNB=tIA->connAtoms.begin();
                    iNB !=tIA->connAtoms.end(); iNB++)
        {
            REAL aOrd = getBondOrder(tMol, tIA->seriNum, *iNB);
            std::cout << "bond order between atom " << tIA->seriNum+1 
                      << " and " << tMol.atoms[*iNB].seriNum+1
                      << " is " << aOrd << std::endl;
            if (aOrd >0)
            {
                tVal +=aOrd;
                std::cout << "total order now " << tVal << std::endl;
            }
            else
            {
                std::cout << "Can not find the bond between atoms " << tIA->id 
                          << " serial number " << tIA->seriNum + 1
                          << " and " << tMol.atoms[*iNB].id 
                          << " serial number " << tMol.atoms[*iNB].seriNum+1
                          << std::endl;
                std::cout << "Some thing is wrong in the Bond list " << std::endl;
                exit(1);
            }
        }
        
        return tVal;
    }
    
    extern REAL getBondOrder(Molecule & tMol, int tIdx1, int tIdx2)
    {
        REAL tOrd = -1.0;
        
        for (std::vector<BondDict>::iterator iB=tMol.bonds.begin();
                iB !=tMol.bonds.end(); iB++)
        {
            if ((iB->atomsIdx[0]==tIdx1 && iB->atomsIdx[1]==tIdx2)
                || (iB->atomsIdx[0]==tIdx2 && iB->atomsIdx[1]==tIdx1))
            {
                tOrd = StrToReal(iB->order);
                if (tOrd ==4.0)
                {
                    tOrd = 1.5;
                }
                
                std::cout << "Bond order " << iB->order << std::endl
                          << tOrd << std::endl;
                break;
            }
        }
        
        return tOrd;
    }
}