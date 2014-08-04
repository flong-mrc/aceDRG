
/* 
 * File:   chiral.cpp
 * Author: flong
 *
 * Created on Oct 4, 2011, 12:29 PM
 */

#include "chiral.h"
#include "atomAssembly.h"

namespace LIBMOL
{
    Chiral::Chiral() :isItTouched(false),
            itsName(NullString),
            itsID(NullString),
            itsSeriNum(ZeroInt),
            itsValue(ZeroReal),
            itsValueSt(ZeroReal),
            itsSigValue(ZeroReal),
            itsValuePre(ZeroReal),
            itsForceConst(ZeroReal)
            
    {
    }
    
    Chiral::Chiral(const Chiral & tC) :isItTouched(tC.isItTouched),
            itsName(tC.getName()),
            itsID(tC.getID()),
            itsSeriNum(tC.getSeriNum()),
            itsValue(tC.getValue(false)),
            itsValueSt(tC.getValue(true)),
            itsSigValue(tC.getSigValue()),
            itsValuePre(tC.getValuePre()),
            itsForceConst(tC.getForceConst())
    {
        for (int i=0; i < (int)tC.atoms.size(); i++)
        {
            atoms.push_back(tC.atoms[i]);
        }
    }
    
    Chiral::~Chiral()
    {
        if(!atoms.empty())
        {
            atoms.clear();
        }
    }

    //Chiral & Chiral::operator=(const Chiral &tC)
    //{
    //}
    
    Name Chiral::getName() const
    {
            return itsName;
    }
    void Chiral::setName(Name tNa)
    {
        itsName = tNa;
    }
        
    ID Chiral::getID() const
    {
        return itsID;
    }
    void Chiral::setID(ID tID)
    {
        itsID = tID;
    }
        
    SeriNumber Chiral::getSeriNum() const
    {
        return itsSeriNum;
    }
    void Chiral::setSeriNum(SeriNumber tSer)
    {
        itsSeriNum = tSer;
    }
    
    REAL Chiral::getValue(bool tL) const
    {
        if(tL)
        {
            return itsValueSt;
        }
        else
        {
            return itsValue;
        }
    }
    void Chiral::setValue(REAL tV, bool tL)
    {
        if(tL)
        {
            itsValueSt = tV;
        }
        else
        {
            itsValue  = tV;
        }
    }
    REAL Chiral::setValue(Atom & tA1, Atom & tA2, Atom & tA3, Atom & tA4)
    {
        /*
        std::cout << " Atom 1 " << tA1.getName() 
                  << " Residue " << tA1.getSeqNum()
                  << " Chain " << tA1.getChainID() 
                  << std::endl;
        std::cout << " Atom 2 " << tA2.getName() 
                  << " Residue " << tA2.getSeqNum()
                  << " Chain " << tA2.getChainID()
                  << std::endl;
        std::cout << " Atom 3 " << tA3.getName() 
                  << " Residue " << tA3.getSeqNum()
                  << " Chain " << tA3.getChainID()
                 << std::endl;
        std::cout << " Atom 4 " << tA4.getName() 
                  << " Residue " << tA4.getSeqNum()
                  << " Chain "   << tA4.getChainID()
                  << std::endl;
        */
        if (tA1.coords.size() ==3 && tA2.coords.size() ==3 
            && tA3.coords.size() ==3 && tA4.coords.size() ==3)
        {
            std::vector<REAL> tV1, tV2, tV3;
            for (int i=0; i < 3; i++)
            {
                tV1.push_back(tA2.coords[i]-tA1.coords[i]);
                tV2.push_back(tA3.coords[i]-tA1.coords[i]);
                tV3.push_back(tA4.coords[i]-tA1.coords[i]);
                
            }
            // std::cout << " getDet " << getDet(tV1, tV2, tV3) << std::endl;
            return getDet(tV1, tV2, tV3);
        }
        else
        {
            return 0.0;
        }
    }
    void Chiral::setValue()
    {
        if (atoms.size() == 4)
        {
            itsValue = setValue(atoms[0], atoms[1], atoms[2], atoms[3]);
        }
        else
        {
            itsValue = 0.0;
        }
    }
    
    REAL Chiral::getSigValue() const
    {
        return itsSigValue;
    }

    void Chiral::setSigValue(REAL tV)
    {
        itsSigValue = tV;
    }
    
    REAL Chiral::getValuePre() const
    {
        return itsValuePre;
    }
    void Chiral::setValuePre(REAL tV)
    {
        itsValuePre = tV;
    }
    
    REAL Chiral::getForceConst() const
    {
        return  itsForceConst;
    }
    void Chiral::setForceConst(REAL tF)
    {
        itsForceConst = tF;
    }
    
    // Another class on Chiral center
    
    ChiralDict::ChiralDict():id(NullString),
            archID(NullString),
            archPos(ZeroInt),
            value(ZeroReal),
            sign(NullString),
            valueST(ZeroReal),
            signST(NullString),
            fConst(100.0)
    {
    }
    
    ChiralDict::ChiralDict(const ChiralDict& tCh):id(tCh.id),
            archID(tCh.archID),
            archPos(tCh.archPos),
            value(tCh.value),
            sign(tCh.sign),
            valueST(tCh.valueST),
            signST(tCh.signST), 
            fConst(tCh.fConst)
    {
        for (std::vector<int>::const_iterator iC=tCh.atoms.begin();
                iC != tCh.atoms.end(); iC++)
        {
            atoms.push_back(*iC);
        }
        
        for (std::map<int, std::vector<int> >::const_iterator iM1=tCh.mutTable.begin();
                iM1 !=tCh.mutTable.end(); iM1++)
        {
            for (std::vector<int>::const_iterator iM2=iM1->second.begin();
                 iM2 != iM1->second.end(); iM2++)
            {
                mutTable[iM1->first].push_back(*iM2);
            }

        }
    }
    
    ChiralDict::~ChiralDict()
    {
    }
    
    void ChiralDict::setMutTable(int tChiralIdx)
    {
        
        
                if((int)atoms.size() ==5)
                { 
                    if (tChiralIdx==1 ||
                        tChiralIdx==2)
                    {
                        mutTable[atoms[1]].push_back(atoms[4]);
                        mutTable[atoms[1]].push_back(atoms[3]);
                        mutTable[atoms[1]].push_back(atoms[2]);

                        mutTable[atoms[2]].push_back(atoms[1]);
                        mutTable[atoms[2]].push_back(atoms[3]);
                        mutTable[atoms[2]].push_back(atoms[4]);

                        mutTable[atoms[3]].push_back(atoms[1]);
                        mutTable[atoms[3]].push_back(atoms[2]);
                        mutTable[atoms[3]].push_back(atoms[4]);

                        mutTable[atoms[4]].push_back(atoms[3]);
                        mutTable[atoms[4]].push_back(atoms[1]);
                        mutTable[atoms[4]].push_back(atoms[2]);
                    }

                    if (tChiralIdx==-1)
                    {
                        mutTable[atoms[1]].push_back(atoms[2]);
                        mutTable[atoms[1]].push_back(atoms[3]);
                        mutTable[atoms[1]].push_back(atoms[4]);

                        mutTable[atoms[2]].push_back(atoms[3]);
                        mutTable[atoms[2]].push_back(atoms[1]);
                        mutTable[atoms[2]].push_back(atoms[4]);

                        mutTable[atoms[3]].push_back(atoms[2]);
                        mutTable[atoms[3]].push_back(atoms[1]);
                        mutTable[atoms[3]].push_back(atoms[4]);

                        mutTable[atoms[4]].push_back(atoms[1]);
                        mutTable[atoms[4]].push_back(atoms[3]);
                        mutTable[atoms[4]].push_back(atoms[2]);
                    }
                }
                else if((int)atoms.size() ==4)
                {
                    if (tChiralIdx==1 ||
                        tChiralIdx==2)
                    {
                        mutTable[atoms[1]].push_back(atoms[3]);
                        mutTable[atoms[1]].push_back(atoms[2]);

                        mutTable[atoms[2]].push_back(atoms[1]);
                        mutTable[atoms[2]].push_back(atoms[3]);

                        mutTable[atoms[3]].push_back(atoms[1]);
                        mutTable[atoms[3]].push_back(atoms[2]);

                    }

                    if (tChiralIdx==-1)
                    {
                        mutTable[atoms[1]].push_back(atoms[2]);
                        mutTable[atoms[1]].push_back(atoms[3]);

                        mutTable[atoms[2]].push_back(atoms[3]);
                        mutTable[atoms[2]].push_back(atoms[1]);

                        mutTable[atoms[3]].push_back(atoms[2]);
                        mutTable[atoms[3]].push_back(atoms[1]);

                    }
                }
                
        
    }
    
    
    void ChiralDict::setMutTable2(int tChiralIdx)
    {
        if((int)atoms.size() ==5)
        { 
            if (tChiralIdx==2 ||
                tChiralIdx==3)
            {
                mutTable[atoms[1]].push_back(atoms[4]);
                mutTable[atoms[1]].push_back(atoms[3]);
                mutTable[atoms[1]].push_back(atoms[2]);

                mutTable[atoms[2]].push_back(atoms[1]);
                mutTable[atoms[2]].push_back(atoms[3]);
                mutTable[atoms[2]].push_back(atoms[4]);

                mutTable[atoms[3]].push_back(atoms[1]);
                mutTable[atoms[3]].push_back(atoms[2]);
                mutTable[atoms[3]].push_back(atoms[4]);

                mutTable[atoms[4]].push_back(atoms[3]);
                mutTable[atoms[4]].push_back(atoms[1]);
                mutTable[atoms[4]].push_back(atoms[2]);
            }
            else if (tChiralIdx==1)
            {
                mutTable[atoms[1]].push_back(atoms[2]);
                mutTable[atoms[1]].push_back(atoms[3]);
                mutTable[atoms[1]].push_back(atoms[4]);

                mutTable[atoms[2]].push_back(atoms[3]);
                mutTable[atoms[2]].push_back(atoms[1]);
                mutTable[atoms[2]].push_back(atoms[4]);

                mutTable[atoms[3]].push_back(atoms[2]);
                mutTable[atoms[3]].push_back(atoms[1]);
                mutTable[atoms[3]].push_back(atoms[4]);

                mutTable[atoms[4]].push_back(atoms[1]);
                mutTable[atoms[4]].push_back(atoms[3]);
                mutTable[atoms[4]].push_back(atoms[2]);
            }
        }
        else if((int)atoms.size() ==4)
        {
            if (tChiralIdx==2 ||
                tChiralIdx==3)
            {
                mutTable[atoms[1]].push_back(atoms[3]);
                mutTable[atoms[1]].push_back(atoms[2]);

                mutTable[atoms[2]].push_back(atoms[1]);
                mutTable[atoms[2]].push_back(atoms[3]);

                mutTable[atoms[3]].push_back(atoms[1]);
                mutTable[atoms[3]].push_back(atoms[2]);
            }
            else if (tChiralIdx==1)
            {
                mutTable[atoms[1]].push_back(atoms[2]);
                mutTable[atoms[1]].push_back(atoms[3]);

                mutTable[atoms[2]].push_back(atoms[3]);
                mutTable[atoms[2]].push_back(atoms[1]);

                mutTable[atoms[3]].push_back(atoms[2]);
                mutTable[atoms[3]].push_back(atoms[1]);

            }
        }       
    }
    
    
    extern int inChirals(std::vector<ChiralDict> tChirals, 
                          AtomDict & tInAtom)
    {
        int tFind =-1;
        for (int iCh=0; iCh <(int)tChirals.size(); iCh++)
        {
            if (tChirals[iCh].archID==tInAtom.id && 
                std::find(tInAtom.inChirals.begin(), tInAtom.inChirals.end(), iCh)
                    ==tInAtom.inChirals.end())
            {
                tInAtom.inChirals.push_back(iCh);
                tFind = iCh;
                break;
            }
        }
        
        return tFind;
    }
    
    extern void buildChiralCluster(std::vector<int>  & tInAts,
                                   std::vector<int>  & tOutAts,
                                   int                 tExc, 
                                   int                 tChStat)
    {
        // check if one of chiral atoms is already in the atom cluster (ring criteria)
        // tInAts are atoms from a chiral set.
        
        // For sp3 
        
        int iPos =-1;
        for (int i=1; i < (int)tInAts.size(); i++)
        {
            if (tInAts[i] !=tExc &&std::find(tOutAts.begin(), tOutAts.end(), tInAts[i])
                    !=tOutAts.end())
            {
                iPos = i;
                break;
            }
        }
        //std::cout << "iPos " << iPos << std::endl;
        
        std::vector<int> tVec;
        for (int j=1; j < (int)tInAts.size(); j++)
        {
            if (j !=iPos && tInAts[j] !=tExc )
            {
                tVec.push_back(tInAts[j]);
                // std::cout << "included " << tInAts[j] << std::endl;
            }
        }
        
       
        
        if (tChStat >0)
        {
            // positive chiral 
            if (iPos==1)
            {
                tOutAts.push_back(tVec[0]);
                if((int)tVec.size()==2)
                {
                    tOutAts.push_back(tVec[1]);
                }
            }
            else if(iPos==2)
            {
                if((int)tVec.size()==2)
                {
                    tOutAts.push_back(tVec[1]);
                }
                tOutAts.push_back(tVec[0]);
            }
            else if(iPos==3)
            {
                tOutAts.push_back(tVec[0]);
                if((int)tVec.size()==2)
                {
                    tOutAts.push_back(tVec[1]);
                }
            }
            else
            {
                tOutAts.push_back(tVec[0]);
                tOutAts.push_back(tVec[1]);
                if((int)tVec.size()==3)
                {
                    tOutAts.push_back(tVec[2]);
                }
            }
        }
        else if (tChStat <0)
        {
            // negative chiral 
            if (iPos==1 || iPos==3)
            {
                if((int)tVec.size()==2)
                {
                    tOutAts.push_back(tVec[1]);
                }
                tOutAts.push_back(tVec[0]);
                
            }
            else if(iPos==2)
            {
                tOutAts.push_back(tVec[0]);
                if((int)tVec.size()==2)
                {
                    tOutAts.push_back(tVec[1]);
                }
            }
            else
            {
                if((int)tVec.size()==3)
                {
                    tOutAts.push_back(tVec[1]);
                    tOutAts.push_back(tVec[2]);
                }
                else
                {
                    tOutAts.push_back(tVec[1]);
                }
                tOutAts.push_back(tVec[0]);
            }
        }
        else
        {
            tOutAts.push_back(tVec[0]);
            tOutAts.push_back(tVec[1]);
            if((int)tVec.size()==3)
            {
                tOutAts.push_back(tVec[2]);
            }
        } 
        /*
        for (int j=0; j < (int)tOutAts.size(); j++)
        {
           
                std::cout << "F included " << tOutAts[j] << std::endl;
            
        }
        */
    }
    
    extern void buildChiralCluster2(ChiralDict        & tChiral,
                                    std::vector<int>  & tOutAts, 
                                    int                 tExc)
    {
        
        
        // For sp3
        
        // check if one of chiral atoms is already in the atom cluster (ring criteria)
        // tInAts are atoms from a chiral set.
        /*
        std::vector<int> tInAts;
        
        // remove the same ring atom 
        for (std::vector<int>::iterator iA=tInAts.begin();
                iA != tInAts.end(); iA++)
        {   
            if (std::find(tOutAts.begin(), tOutAts.end(), *iA)!=tOutAts.end())
            { 
                tInAts.erase(iA);
                std::cout << (int)tInAts.size() << std::endl;
                break;
            }
        }
        
        
        // remove one bonding atom
        std::vector<int>::iterator iFind 
                   = std::find(tInAts.begin(), tInAts.end(), tExc);
        if (iFind!=tInAts.end())
        {
             tInAts.erase(iFind);
        }
        
        // Now atom list associated with this chiral center contain
        // only atoms need to put into the output atom list
        
        if (tChStat >0)
        {
            // positive chiral
            if ((int)tInAts.size() ==3)
            {
                // no deleted atom from the input list, use the original sequence.
                for (std::vector<int>::iterator iAt=tInAts.begin();
                        iAt != tInAts.end(); iAt++)
                {
                    tOutAts.push_back(*iAt);
                }
            }
            else
            {
                // some atoms are deleted, using reverse sequence 
                // of  the input list
                for (std::vector<int>::reverse_iterator iAt=tInAts.rbegin();
                        iAt !=tInAts.rend(); iAt++)
                {
                    tOutAts.push_back(*iAt);
                }
            }
        }
        else
        {
            // negative chiral center
            if ((int)tInAts.size() ==3)
            {
                // no deleted atom from the input list, use the reverse sequence.
                for (std::vector<int>::reverse_iterator iAt=tInAts.rbegin();
                        iAt != tInAts.rend(); iAt++)
                {
                    tOutAts.push_back(*iAt);
                }
            }
            else
            {
                // some atoms are deleted, using original sequence 
                // of  the input list
                for (std::vector<int>::iterator iAt=tInAts.begin();
                        iAt !=tInAts.end(); iAt++)
                {
                    tOutAts.push_back(*iAt);
                }
            }
            
        }
        */
        if((int)tOutAts.size()==0)
        {
            // std::cout << "Here " << std::endl;
            
            for (std::vector<int>::iterator iA = tChiral.mutTable[tExc].begin();
                 iA != tChiral.mutTable[tExc].end(); iA++)
            {
                if (std::find(tOutAts.begin(), tOutAts.end(), *iA)==tOutAts.end())
                {
                    tOutAts.push_back(*iA);
                    // std::cout << *iA << std::endl;        
                }
            }
        }
        else
        {
            getMutList(tChiral.mutTable[tExc], tOutAts);
        }
    }
    
    extern void getMutList(std::vector<int> & tIn, std::vector<int> & tOut)
    {
        int tPos =-1;
        if ((int)tOut.size()==1)
        {
            for (int i=0; i < (int)tIn.size(); i++)
            {
                if (tIn[i]==tOut[0])
                {
                    tPos = i;
                    break;
                }
            }
        }
        
       
        if (tPos==-1)
        {
            for (std::vector<int>::iterator iP=tIn.begin();
                        iP !=tIn.end(); iP++)
            {
                tOut.push_back(*iP);
                //std::cout << "tOut " << *iP << std::endl;
            }
        }
        else
        {
            //std::cout << "tPOs " << tPos << std::endl;
            int tSize = (int)tIn.size()-1;
            for(int i=0; i < tSize; i++)
            {
                int j = tPos+1+i;
                if (j > tSize)
                {
                    j=j-(int)tIn.size();
                }
                tOut.push_back(tIn[j]);
                // std::cout << "tIn " << j << " = " << tIn[j] << std::endl;
            }
        }
    }
}