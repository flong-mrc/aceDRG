
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
            if (inFile.is_open())
            {
                setupSystem();
            }
            else
            {
                std::cout << tFname << " can not be open for reading. Check the file "
                         << std::endl;
                exit(1);
            }
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
              containMetal(false),
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
                                // std::cout << "aAtom.seriNum " << tAtomIdx << std::endl;
                                aAtom.coords[0] = StrToReal(tBuf[0]);
                                aAtom.coords[1] = StrToReal(tBuf[1]);
                                aAtom.coords[2] = StrToReal(tBuf[2]);
                                aAtom.coordExist = true;
                                aAtom.chemType  = TrimSpaces(tBuf[3]);
                                aAtom.formalCharge = strToCharge(tBuf[5]);
                                aAtom.inChiralIdx = StrToInt(tBuf[6]); 
                                aAtom.chiralIdx = StrToInt(tBuf[6]); 
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
                                aBond.order = tBuf[2];
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
            setHAtomCoordsMols(i);
           
            
            for (std::vector<BondDict>::iterator iB=allMols[i].bonds.begin();
                    iB !=allMols[i].bonds.end(); iB++)
            {
                std::string tS = iB->order;
                OrderStrToStr(tS, iB->order);
            }
            allMols[i].hasCoords = false;
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
                          << " charge " << iA->formalCharge << std::endl
                          << "Coordinates : " 
                          << iA->coords[0] << "\t" 
                          << iA->coords[1] << "\t"
                          << iA->coords[2] << std::endl;
                if (iA->connAtoms.size() !=0)
                {
                    std::cout << "It connects to the following atoms: " <<std::endl;
                    for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                             iNB !=iA->connAtoms.end(); iNB++)
                    {
                        std::cout << "Atom " << iM->atoms[*iNB].id << std::endl;
                    }
                }
                if (iA->connHAtoms.size() !=0)
                {
                    std::cout << "It bonds to " << iA->connHAtoms.size()  
                              << "H atoms " << std::endl;
                }
                std::cout << "Its is in " << iA->inChirals.size() << " chirals " << std::endl;
                if (iA->inChirals.size())
                {
                    
                    for (std::vector<int>::iterator iCh=iA->inChirals.begin();
                            iCh!=iA->inChirals.end(); iCh++)
                    {
                        std::cout << "For chiral center " << *iCh << std::endl;
                        std::cout << "Its mutational table is " << std::endl;
                        for (std::map<int, std::vector<int> >::iterator iMu=iM->chirals[*iCh].mutTable.begin(); 
                             iMu !=iM->chirals[*iCh].mutTable.end(); iMu++)
                        {
                            std::cout << "start from atom " << iM->atoms[iMu->first].id 
                                      << " and " << iA->id << " mutable is " << std::endl;
                            for (std::vector<int>::iterator iMu2=iMu->second.begin();
                                   iMu2 != iMu->second.end(); iMu2++)
                            {
                                std::cout << iM->atoms[*iMu2].id << "\t";
                            }
                            std::cout << std::endl;
                        }
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
        
        std::vector<std::string> orgElems;
        orgElems.push_back("C");
        orgElems.push_back("O");
        orgElems.push_back("N");
        orgElems.push_back("S");
        orgElems.push_back("P");
        orgElems.push_back("B");
        
        orgElems.push_back("Si");
        orgElems.push_back("S");
        
        PeriodicTable aPTab;
        
        for (std::vector<AtomDict>::iterator iA=allMols[tIdxMol].atoms.begin();
                iA !=allMols[tIdxMol].atoms.end(); iA++)
        {
            
            if (std::find(orgElems.begin(), orgElems.end(), iA->chemType)
                 !=orgElems.end())
            {
                // need to check protonated form
                REAL addH =checkProtonateAll(iA, allMols[tIdxMol], aPTab);
                // REAL addH =checkProtonated(iA, tIdxMol);
                if (addH >0)
                {
                    addHAtoms(tIdxMol, iA->seriNum, addH);
                }
            }
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
        
        reNameHAtoms(tIdxMol);
    }
    
    void MolSdfFile::setHAtomCoordsMols(int tIdxMol)
    {
        
        for (std::vector<AtomDict>::iterator iA=allMols[tIdxMol].atoms.begin();
                iA !=allMols[tIdxMol].atoms.end(); iA++)
        {
            if (iA->chemType=="H" && iA->coordExist==false)
            {
                if (iA->connAtoms.size() ==1 )
                {
                    if (allMols[tIdxMol].atoms[iA->connAtoms[0]].bondingIdx==3)
                    {
                        
                        setOneHAtomCoordsSP3(allMols[tIdxMol].atoms, iA);
                    }
                    else if (allMols[tIdxMol].atoms[iA->connAtoms[0]].bondingIdx==2)
                    {
                        setOneHAtomCoordsSP2(allMols[tIdxMol].atoms, iA);
                    }
                    else if (allMols[tIdxMol].atoms[iA->connAtoms[0]].bondingIdx==1)
                    {
                        setOneHAtomCoordsSP(allMols[tIdxMol].atoms, iA);
                    }
                    else
                    {
                        std::cout << "The bondingIdx for atom " 
                                  <<  allMols[tIdxMol].atoms[iA->connAtoms[0]].id
                                  << " is " << allMols[tIdxMol].atoms[iA->connAtoms[0]].bondingIdx
                                  << std::endl;
                    }
                }
                else
                {
                    std::cout << "Bug. H atom " << iA->id << " connects to "
                              << iA->connAtoms.size() << " atom(s)" << std::endl;
                    exit(1);
                }
            }
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
            
            std::string tS;
            getDigitSec(allMols[tIdxMol].atoms[tIdxAtm].id, tS);
            if (allMols[tIdxMol].atoms[tIdxAtm].chemType=="C")
            {
                aH.id  = aH.chemType + tS +IntToStr(i+1);
            }
            else
            {
                aH.id  = aH.chemType + allMols[tIdxMol].atoms[tIdxAtm].id;
            }
            
            
            aH.seriNum =  (int)allMols[tIdxMol].atoms.size()
                          +(int)allMols[tIdxMol].extraHAtoms.size();
            aH.connAtoms.push_back(allMols[tIdxMol].atoms[tIdxAtm].seriNum);
            allMols[tIdxMol].atoms[tIdxAtm].connAtoms.push_back(aH.seriNum);
            allMols[tIdxMol].atoms[tIdxAtm].connHAtoms.push_back(aH.seriNum);
            // allMols[tIdxMol].atoms.push_back(aH);
            allMols[tIdxMol].extraHAtoms.push_back(aH);
            
        }
        // allMols[tIdxMol].hasCoords = false;
        
        std::cout << tNumH << " H atoms have been added to bind atom " 
                  << allMols[tIdxMol].atoms[tIdxAtm].id << std::endl;  
    }
    
    void MolSdfFile::reNameHAtoms(int tIdxMol)
    {
        int aIdx =1;
        for (std::vector<AtomDict>::iterator iHA=allMols[tIdxMol].atoms.begin();
                iHA !=allMols[tIdxMol].atoms.end(); iHA++)
        {
            if (iHA->chemType.compare("H")==0)
            {
                iHA->id = "H" + IntToStr(aIdx);
                aIdx+=1;
            }
        }
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
                    // iAt->chiralIdx  = 1;
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
                    if (iAt->chiralIdx==0)
                    {
                        iAt->chiralIdx  = 1;
                    }
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
                else if (t_len==2)
                {
                    iAt->bondingIdx = 2;
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
                // iA->inChirals.push_back(tChIdx);
                ChiralDict aCh;
                
                aCh.atoms.push_back(iA->seriNum);
               
                // aCh.seriNum = (int)allMols[tIdxMol].chirals.size();
                
                
                
                //std::cout << "Atom is " << iA->id << std::endl;
                //std::cout << "Its serial number " << iA->seriNum << std::endl;
                //std::cout << "number of bonds " << allMols[tIdxMol].bonds.size() 
                //          << std::endl;
                 
                // int i=0;
                std::vector<int> hIdx;
                for (std::vector<BondDict>::iterator iB=allMols[tIdxMol].bonds.begin();
                        iB !=allMols[tIdxMol].bonds.end(); iB++)
                {
                    //if (i <3)
                    //{
                        int idx0=iB->atomsIdx[0], idx1=iB->atomsIdx[1];
                        // std::cout << "Bond atoms: " << allMols[tIdxMol].atoms[idx0].seriNum 
                        //          << " and " << allMols[tIdxMol].atoms[idx1].seriNum  << std::endl;
                        if (std::find(iB->atomsIdx.begin(), iB->atomsIdx.end(), iA->seriNum) !=iB->atomsIdx.end() )
                        {
                            //std::cout << "Atom " << iA->id << " is in this bond " << std::endl;   
                            if (idx0 ==iA->seriNum)
                            {
                                if(allMols[tIdxMol].atoms[idx1].chemType !="H")
                                {
                                    aCh.atoms.push_back(idx1);
                                    // allMols[tIdxMol].atoms[idx1].inChirals.push_back(tChIdx);
                                    //std::cout << allMols[tIdxMol].atoms[idx1].id << " is included 1 " << std::endl;
                                    //i++;
                                }
                                else
                                {
                                    hIdx.push_back(idx1);
                                }
                            }
                            else if (idx1 ==iA->seriNum)
                            {
                                if(allMols[tIdxMol].atoms[idx0].chemType !="H")
                                {
                                    aCh.atoms.push_back(idx0);
                                    //std::cout << allMols[tIdxMol].atoms[idx0].id << " is included 0 " << std::endl;
                                    //i++;
                                }
                                else
                                {
                                    hIdx.push_back(idx1);
                                }
                            }
                        }
                    //}
                }
                if (hIdx.size() > 0)
                {
                    for (std::vector<int>::iterator iH=hIdx.begin();
                            iH !=hIdx.end(); iH++)
                    {
                        aCh.atoms.push_back(*iH);
                    }
                }
                /*
                for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                        iNB !=iA->connAtoms.end(); iNB++)
                {
                    if (i < 3 && allMols[tIdxMol].atoms[*iNB].chemType !="H")
                    {
                        aCh.atoms.push_back(*iNB);
                        allMols[tIdxMol].atoms[*iNB].inChirals.push_back(tChIdx);
                        i++;
                    }
                }
                 */
                
                //std::cout << "Number of atoms in the chiral center " << aCh.atoms.size() << std::endl;
                
                if (aCh.atoms.size() > 3)
                {
                    // ChiToStr(iA->inChiralIdx, aCh.sign);
                    aCh.sign = "both";
                    aCh.setMutTable2(iA->inChiralIdx);
                    int aNum = (int)allMols[tIdxMol].chirals.size() + 1;
                    aCh.id = "chi" + IntToStr(aNum);
                    
                    iA->inChirals.push_back(tChIdx);
                    
                    allMols[tIdxMol].chirals.push_back(aCh);
                    //std::cout << "Atom " << iA->id << " is the chiral center " 
                    //          << iA->inChirals[0]
                    //          << std::endl;
                }
                
            }
        }
        
    }
    
    SYBLMol2File::SYBLMol2File():hasCoords(true),
            hasConnect(false),
            hasH(false),
            containMetal(false)
    {
    }
    
    SYBLMol2File::SYBLMol2File(Name tFname, std::ios_base::openmode tOpenMode) :
              hasCoords(true),
              hasConnect(false),
              hasH(false),
              containMetal(false)
    {
        if (tOpenMode == std::ios::in)
        {
            inFile.open(tFname.c_str(), tOpenMode);
            if (inFile.is_open())
            {
                setupSystem();
            }
            else
            {
                std::cout << tFname << " can not be open for reading. Check the file "
                         << std::endl;
                exit(1);
            }
        }
        else
        {
            outFile.open(tFname.c_str(), tOpenMode);
        }    
    }
    
    SYBLMol2File::SYBLMol2File(FileName tFname, std::ios_base::openmode tOpenMode) :
              hasCoords(false),
              hasConnect(false),
              hasH(false),
              containMetal(false)
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
    
    SYBLMol2File::~SYBLMol2File()
    {
    }
    
    void SYBLMol2File::setupSystem()
    {
        if (inFile.is_open() )
        {  
            // related to records
             
            std::string tRecord="";
            
            while(!inFile.eof())
            {
                std::getline(inFile, tRecord);
                tRecord = TrimSpaces(tRecord);
                
                if (tRecord.size() !=0)
                {
                    if (tRecord.find("@") !=std::string::npos)
                    {
                        setBlock(tRecord);
                    }
                    else
                    {
                        if (mol2Dict["ATOM"])
                        {
                            getAtomInfo(tRecord);
                        }
                        else if (mol2Dict["BOND"])
                        {
                            getBondInfo(tRecord);
                        }
                        //else if(mol2Dict["CRYSIN"])
                        //{
                        //    getCrysInfo(tRecord);
                        //}
                    }
                }
            }
            inFile.close();
            
            
            addAllHAtoms();
            setAtomsBondingAndChiralCenter();
            setChiral();
            
            
            
            if (atoms.size() !=0)
            {
                std::cout << "The system contains " << atoms.size() << std::endl
                      << "Including " << extraHAtoms.size() << " H atoms " 
                      << std::endl;
                
                for (std::vector<AtomDict>::iterator iAt=atoms.begin();
                        iAt != atoms.end(); iAt++)
                {
                    std::cout << std::endl << "For atom " << iAt->seriNum
                              << ": its name is " << iAt->id
                              << ", its element symbol is "
                              << iAt->chemType << std::endl;
                    if (iAt->connAtoms.size() >0)
                    {
                        std::cout << "It bonds to " << iAt->connAtoms.size()
                                  << " atoms. These atoms are: " << std::endl;
                        for (std::vector<int>::iterator iNB=iAt->connAtoms.begin();
                                iNB !=iAt->connAtoms.end(); iNB++)
                        {
                            std::cout << "Atom " << atoms[*iNB].id << std::endl;
                        }
                    }
                }
            }
            else
            {
                std::cout << " Bug or wrong file format. system contain no atoms "
                          << std::endl;
            }
            
            if (bonds.size() !=0)
            {
                std::cout << "The system contains " << bonds.size() 
                          << " bonds. They are: " <<  std::endl;
                
                for (std::vector<BondDict>::iterator iBo=bonds.begin();
                        iBo !=bonds.end(); iBo++)
                {
                    std::cout << "Bond " << iBo->seriNum;
                    
                    if (iBo->atomsIdx.size()==2)
                    {
                      std::cout << ", which consists of atom, "
                                << iBo->atoms[0] << " and "
                                << iBo->atoms[1] << std::endl;
                    }
                    else
                    {
                        std::cout << "Bug, the bond consists of less than 2 atoms "
                                  << std::endl;
                    }
                }
            }
            else
            {
                std::cout << "Warning: check the input mol2 file. it contains no bonds"
                          << std::endl;
            } 
            
        }
        
        
    }
    
    void SYBLMol2File::setBlock(std::string tLine)
    {
        StrUpper(tLine);
        if (tLine.find("ATOM") !=std::string::npos )
        {
            mol2Dict["ATOM"] = true;
            for (std::map<ID, bool>::iterator iM=mol2Dict.begin();
                    iM !=mol2Dict.end(); iM++)
            {
                if (iM->first !="ATOM")
                {
                    iM->second = false; 
                }
            }
        }
        else if(tLine.find("BOND") !=std::string::npos)
        {
            mol2Dict["BOND"] = true;
            for (std::map<ID, bool>::iterator iM=mol2Dict.begin();
                    iM !=mol2Dict.end(); iM++)
            {
                if (iM->first !="BOND")
                {
                    iM->second = false; 
                } 
            }
        }
        else if(tLine.find("CRYSIN") !=std::string::npos)
        {
            mol2Dict["CRYSIN"] = true;
            for (std::map<ID, bool>::iterator iM=mol2Dict.begin();
                    iM !=mol2Dict.end(); iM++)
            {
                if (iM->first !="CRYSIN")
                {
                    iM->second = false; 
                }
            }
        }
        else
        {
            for (std::map<ID, bool>::iterator iM=mol2Dict.begin();
                    iM !=mol2Dict.end(); iM++)
            {
                iM->second = false; 
                
            }
        }
    }
  
    
    void SYBLMol2File::iniDict()
    {
        mol2Dict["ATOM"]     = false;
        mol2Dict["BOND"]     = false;
        mol2Dict["ALT_TYPE"] = false;
        mol2Dict["CRYSIN"]   = false;
    }
    
    void SYBLMol2File::getAtomInfo(std::string tLine)
    {
        std::vector<std::string> tStrs;
        StrTokenize(tLine, tStrs);
        if (tStrs.size()>=6)
        {
            AtomDict  aAtom;
            aAtom.seriNum = StrToInt(tStrs[0])-1;
            aAtom.id      = tStrs[1];
            aAtom.coords.push_back(StrToReal(tStrs[2]));
            aAtom.coords.push_back(StrToReal(tStrs[3]));
            aAtom.coords.push_back(StrToReal(tStrs[4]));
            std::vector<std::string> tStrs1;
            StrTokenize(tStrs[5], tStrs1, '.');
            if (tStrs1.size() >0)
            {
                aAtom.chemType = tStrs1[0];
                if (tStrs1.size() > 1)
                {
                    atomSYBYLTypes[aAtom.id] = tStrs1[1];
                }
            }
            else
            {
                std::cout << "ATOM line format error: atom type field is empty "
                          << std::endl;
                exit(1);
            }
            if (tStrs.size()>=8)
            {
                aAtom.resName = tStrs[7];
                if (tStrs.size()>=9)
                {
                    aAtom.charge = StrToReal(tStrs[8]);
                }
            }
            
            atoms.push_back(aAtom);
            std::cout << "Atom line "<< tLine << std::endl;
            std::cout << " Set an atom: seriNum : " << aAtom.seriNum 
                      << ", id : " << aAtom.id << " element symbol : " 
                      << aAtom.chemType;
            
            if (aAtom.resName.size() !=0)
            {
                std::cout << " residueID : " << aAtom.resName;
            }
            if (atomSYBYLTypes.find(aAtom.id) !=atomSYBYLTypes.end())
            {
                std::cout << " atom SYBYL type : " << atomSYBYLTypes[aAtom.id]; 
            }
            std::cout << std::endl;
        }
    }
   
    void SYBLMol2File::getBondInfo(std::string tLine)
    {
        std::vector<std::string> tStrs;
        StrTokenize(tLine, tStrs);
        
        if (tStrs.size() >=4)
        {
            BondDict aBond;
            aBond.seriNum = StrToInt(tStrs[0])-1;
            aBond.atomsIdx.push_back(StrToInt(tStrs[1])-1);
            aBond.atomsIdx.push_back(StrToInt(tStrs[2])-1);
            aBond.order   = TrimSpaces(tStrs[3]);
            aBond.orderN  = StrToOrder2(aBond.order);
            
            int i=0;
            for (std::vector<AtomDict>::iterator iAt=atoms.begin();
                    iAt !=atoms.end(); iAt++)
            {
                if (iAt->seriNum==aBond.atomsIdx[0])
                {
                    iAt->connAtoms.push_back(aBond.atomsIdx[1]);
                    aBond.atoms.push_back(iAt->id);
                    aBond.fullAtoms[iAt->id] = iAt->seriNum;
                    i++;
                }
                else if (iAt->seriNum==aBond.atomsIdx[1])
                {
                    iAt->connAtoms.push_back(aBond.atomsIdx[0]);
                    aBond.atoms.push_back(iAt->id);
                    aBond.fullAtoms[iAt->id] = iAt->seriNum;
                    i++;
                }
                
                if (i==2)
                {
                    break;
                }
            }
            bonds.push_back(aBond);
            // check 
            std::cout << "A bond line " << tLine << std::endl;
            std::cout << "A bond is added. Its seriNum " 
                      << aBond.seriNum << ", contain "
                      << aBond.atomsIdx.size() << " atoms " 
                      << " bond order is " << aBond.order 
                      << " bond order value is " << aBond.orderN 
                      << std::endl;
                      
            if (aBond.atomsIdx.size() == 2)
            {
                std::cout << "Two atoms are: atom " << aBond.atomsIdx[0]
                          << " id " << atoms[aBond.atomsIdx[0]].id
                          << " and atom " << aBond.atomsIdx[1]
                          << " id " << atoms[aBond.atomsIdx[1]].id 
                          << std::endl;
            }
            else 
            {
                std::cout << "Bug. The bond contains less than 2 atom " 
                          << std::endl;  
            }
            
        }
        
    }
    
    void SYBLMol2File::getCrysInfo(std::string tLine)
    {
        // std::vector<std::string> tStrs;
        // StrTokenize(tLine, tStrs);
    }
    
    void SYBLMol2File::addAllHAtoms()
    {
                
        std::vector<std::string> orgElems;
        orgElems.push_back("C");
        orgElems.push_back("O");
        orgElems.push_back("N");
        orgElems.push_back("S");
        orgElems.push_back("P");
        orgElems.push_back("B");
        
        orgElems.push_back("Si");
        orgElems.push_back("S");
        
        PeriodicTable aPTab;
        
        for (std::vector<AtomDict>::iterator iA=atoms.begin();
                iA !=atoms.end(); iA++)
        {
            
            if (std::find(orgElems.begin(), orgElems.end(), iA->chemType)
                 !=orgElems.end())
            {
                // need to check protonated form
                REAL addH =checkProtonateAll(iA, atoms, bonds, aPTab);
                // REAL addH =checkProtonated(iA, tIdxMol);
                if (addH >0)
                {
                    addHAtoms(iA->seriNum, addH);
                }
            }
        }
        
        // Now add newly generated H atoms to atoms in the molecule
        // and add associated bonds to the bond set in the molecule
        for (std::vector<AtomDict>::iterator iHA=extraHAtoms.begin();
                iHA !=extraHAtoms.end(); iHA++)
        {
            atoms.push_back(*iHA);
            BondDict aBo;
            aBo.atomsIdx.push_back(iHA->connAtoms[0]);
            aBo.atomsIdx.push_back(iHA->seriNum);
            aBo.order = "1";
            aBo.seriNum = (int)bonds.size();
            aBo.atoms.push_back(atoms[iHA->connAtoms[0]].id);
            aBo.atoms.push_back(iHA->id);
            aBo.fullAtoms[iHA->id] = iHA->seriNum;
            aBo.fullAtoms[atoms[iHA->connAtoms[0]].id] = atoms[iHA->connAtoms[0]].seriNum;
                                                         bonds.push_back(aBo);
        }
        
        // Rename H atoms
        int aIdx =1;
        for (std::vector<AtomDict>::iterator iHA=atoms.begin();
                iHA !=atoms.end(); iHA++)
        {
            if (iHA->chemType.compare("H")==0)
            {
                iHA->id = "H" + IntToStr(aIdx);
                aIdx+=1;
            }
        } 
        
    }
    
    void SYBLMol2File::addHAtoms(int tIdxAtm, REAL tNumH)
    {
        for (int i=0; i < (int)tNumH; i++)
        {
            AtomDict aH;
            aH.chemType = "H";
            
            std::string tS;
            getDigitSec(atoms[tIdxAtm].id, tS);
            if (atoms[tIdxAtm].chemType=="C")
            {
                aH.id  = aH.chemType + tS +IntToStr(i+1);
            }
            else
            {
                aH.id  = aH.chemType + atoms[tIdxAtm].id;
            }
            
            
            aH.seriNum =  (int)atoms.size()
                         +(int)extraHAtoms.size();
            aH.connAtoms.push_back(atoms[tIdxAtm].seriNum);
            atoms[tIdxAtm].connAtoms.push_back(aH.seriNum);
            atoms[tIdxAtm].connHAtoms.push_back(aH.seriNum);
            // allMols[tIdxMol].atoms.push_back(aH);
            extraHAtoms.push_back(aH);
            
        }
        // allMols[tIdxMol].hasCoords = false;
        
        std::cout << tNumH << " H atoms have been added to bind atom " 
                  << atoms[tIdxAtm].id << std::endl;  
    }
    
    int SYBLMol2File::getNumOxyConnect(std::vector<AtomDict>::iterator iA)
    {
        int nO=0;
        for (std::vector<int>::iterator iC=iA->connAtoms.begin();
                iC !=iA->connAtoms.end(); iC++)
        {
            if (atoms[*iC].chemType.compare("O")==0)
            {
                nO++;
            }
        }
        return nO;
    }
    
    void SYBLMol2File::setAtomsBondingAndChiralCenter()
    {
        // First round
        for (std::vector<AtomDict>::iterator iAt = atoms.begin();
                iAt != atoms.end(); iAt++)
        {
            int t_len =0;
            for (std::vector<int>::iterator iConn=iAt->connAtoms.begin();
                    iConn !=iAt->connAtoms.end(); iConn++)
            {
                if(!atoms[*iConn].isMetal)
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
                    if (getNumOxyConnect(iAt)==1)
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
                    // iAt->chiralIdx  = 1;
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
                    if (iAt->chiralIdx==0)
                    {
                        iAt->chiralIdx  = 1;
                    }
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
                            if (std::find(atps.begin(), atps.end(), atoms[*iNA].chemType)==atps.end())
                            {
                                atps.push_back(atoms[*iNA].chemType);
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
                else if (t_len==2)
                {
                    iAt->bondingIdx = 2;
                }
            }
        }
        
        // more conditions 
        
        for (std::vector<AtomDict>::iterator iAt = atoms.begin();
                iAt != atoms.end(); iAt++)
        {   
            int t_len =0;
            for (std::vector<int>::iterator iConn=iAt->connAtoms.begin();
                    iConn !=iAt->connAtoms.end(); iConn++)
            {
                if(!atoms[*iConn].isMetal)
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
                        if(atoms[*iCA].bondingIdx == 2)
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
        for (std::vector<AtomDict>::iterator iA=atoms.begin();
                iA != atoms.end(); iA++)
        {
            if (iA->chiralIdx !=0)
            {
                
                std::vector<ID> chirRAtms;
                for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                        iNB != iA->connAtoms.end(); iNB++)
                {
                    std::size_t tFind = atoms[*iNB].chemType.find("H");
                    if (tFind !=std::string::npos)
                    {
                        chirRAtms.push_back(atoms[*iNB].id);
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
    
    void SYBLMol2File::setChiral()
    {
        for (std::vector<AtomDict>::iterator iA=atoms.begin();
             iA !=atoms.end(); iA++)
        {
            if (iA->chiralIdx !=0)
            {
                int tChIdx = (int)chirals.size();
                // iA->inChirals.push_back(tChIdx);
                ChiralDict aCh;
                
                aCh.atoms.push_back(iA->seriNum);
               
                // aCh.seriNum = (int)allMols[tIdxMol].chirals.size();
                
                
                
                //std::cout << "Atom is " << iA->id << std::endl;
                //std::cout << "Its serial number " << iA->seriNum << std::endl;
                //std::cout << "number of bonds " << allMols[tIdxMol].bonds.size() 
                //          << std::endl;
                 
                // int i=0;
                std::vector<int> hIdx;
                for (std::vector<BondDict>::iterator iB=bonds.begin();
                        iB !=bonds.end(); iB++)
                {
                    //if (i <3)
                    //{
                        int idx0=iB->atomsIdx[0], idx1=iB->atomsIdx[1];
                        // std::cout << "Bond atoms: " << allMols[tIdxMol].atoms[idx0].seriNum 
                        //          << " and " << allMols[tIdxMol].atoms[idx1].seriNum  << std::endl;
                        if (std::find(iB->atomsIdx.begin(), iB->atomsIdx.end(), iA->seriNum) !=iB->atomsIdx.end() )
                        {
                            //std::cout << "Atom " << iA->id << " is in this bond " << std::endl;   
                            if (idx0 ==iA->seriNum)
                            {
                                if(atoms[idx1].chemType !="H")
                                {
                                    aCh.atoms.push_back(idx1);
                                    // allMols[tIdxMol].atoms[idx1].inChirals.push_back(tChIdx);
                                    //std::cout << allMols[tIdxMol].atoms[idx1].id << " is included 1 " << std::endl;
                                    //i++;
                                }
                                else
                                {
                                    hIdx.push_back(idx1);
                                }
                            }
                            else if (idx1 ==iA->seriNum)
                            {
                                if(atoms[idx0].chemType !="H")
                                {
                                    aCh.atoms.push_back(idx0);
                                    //std::cout << allMols[tIdxMol].atoms[idx0].id << " is included 0 " << std::endl;
                                    //i++;
                                }
                                else
                                {
                                    hIdx.push_back(idx1);
                                }
                            }
                        }
                    //}
                }
                if (hIdx.size() > 0)
                {
                    for (std::vector<int>::iterator iH=hIdx.begin();
                            iH !=hIdx.end(); iH++)
                    {
                        aCh.atoms.push_back(*iH);
                    }
                }
                /*
                for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                        iNB !=iA->connAtoms.end(); iNB++)
                {
                    if (i < 3 && allMols[tIdxMol].atoms[*iNB].chemType !="H")
                    {
                        aCh.atoms.push_back(*iNB);
                        allMols[tIdxMol].atoms[*iNB].inChirals.push_back(tChIdx);
                        i++;
                    }
                }
                 */
                
                //std::cout << "Number of atoms in the chiral center " << aCh.atoms.size() << std::endl;
                
                if (aCh.atoms.size() > 3)
                {
                    // ChiToStr(iA->inChiralIdx, aCh.sign);
                    aCh.sign = "both";
                    aCh.setMutTable2(iA->inChiralIdx);
                    int aNum = (int)chirals.size() + 1;
                    aCh.id = "chi" + IntToStr(aNum);
                    
                    iA->inChirals.push_back(tChIdx);
                    
                    chirals.push_back(aCh);
                    //std::cout << "Atom " << iA->id << " is the chiral center " 
                    //          << iA->inChirals[0]
                    //          << std::endl;
                }
                
            }
        }
    }
    
    
}