/* 
 * File:   Torsion.cpp
 * Author: flong
 *
 * Created on September 12, 2012, 11:32 AM
 */

#include "torsion.h"

namespace LIBMOL
{
    
    class Atom;
    class AtomDict;
    class Bond;
    class BondDict;
    class Angle;
    
    Torsion::Torsion(): varFlag(ZeroInt),
            isItTouched(false),
            itsName(NullString),
            itsID(NullString),
            itsSeriNum(ZeroInt),
            itsValue(ZeroReal),
            itsValueSt(ZeroReal),
            itsSigValue(ZeroReal),
            itsValuePre(ZeroReal),
            itsForceConst(ZeroReal),
            itsPeriod(ZeroInt)
           
    {
    }
    
    Torsion::Torsion(const Torsion & tT):varFlag(tT.varFlag),
            isItTouched(tT.isItTouched),
            itsName(tT.getName()),
            itsID(tT.getID()),
            itsSeriNum(tT.getSeriNum()),
            itsValue(tT.getValue(false)),
            itsValueSt(tT.getValue(true)),
            itsSigValue(tT.getSigValue()),
            itsValuePre(tT.getValuePre()),
            itsForceConst(tT.getForceConst()),
            itsPeriod(tT.getPeriod())
    {   
        for (std::vector<REAL>::const_iterator iR =tT.normalUnit.begin();
                iR != tT.normalUnit.end(); iR++)
        {
            normalUnit.push_back((*iR));
        }
        
        for (std::vector<int>::const_iterator iN =tT.torsUpAtomList.begin();
                iN != tT.torsUpAtomList.end(); iN++)
        {
            torsUpAtomList.push_back((*iN));
        }
        
        for (std::vector<Atom>::const_iterator iA =tT.atoms.begin();
                iA != tT.atoms.end(); iA++)
        {
            atoms.push_back((*iA));
        }
    }
 
  
    Torsion::~Torsion()
    {
        if(!normalUnit.empty())
        {
            normalUnit.clear();
        }
        if(!torsUpAtomList.empty())
        {
            torsUpAtomList.clear();
        }
        if(!atoms.empty())
        {
            atoms.clear();
        }
    }
    
    Name Torsion::getName() const
    {
            return itsName;
    }
    void Torsion::setName(Name tNa)
    {
        itsName = tNa;
    }
        
    ID Torsion::getID() const
    {
        return itsID;
    }
    void Torsion::setID(ID tID)
    {
        itsID = tID;
    }
        
    SeriNumber Torsion::getSeriNum() const
    {
        return itsSeriNum;
    }
    void Torsion::setSeriNum(SeriNumber tSer)
    {
        itsSeriNum = tSer;
    }
           
    REAL Torsion::getValue(bool tB) const
    {
        if(tB)
        {
            return itsValueSt;
        }
        
        return itsValue;
    }
    void Torsion::setValue(REAL tT, bool tB)
    {
        if (tB)
        {
            itsValueSt = tT;
        }
        else 
        {
            itsValue   = tT;
        }
    }
    REAL Torsion::setValue(Atom & tA1, Atom & tA2, Atom & tA3, Atom & tA4)
    {
        std::vector<REAL> tR1, tR2, tR3;
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
        if(tA1.coords.size() == 3 && tA2.coords.size() == 3 && 
                tA3.coords.size() ==3)
        {
            for (int i = 0; i < (int)tA1.coords.size(); i++)
            {
                tR1.push_back(tA2.coords[i]-tA1.coords[i]);
                //std::cout << "tR1 " << i+1 << "  "
                //        <<  tA2.coords[i]-tA1.coords[i] << std::endl;
                tR2.push_back(tA3.coords[i]-tA2.coords[i]);
                //std::cout << "tR2 " << i+1 << "  "
                //        <<  tA3.coords[i]-tA2.coords[i] << std::endl;
                tR3.push_back(tA4.coords[i]-tA3.coords[i]);
                //std::cout << "tR3 " << i+1 << "  "
                //        <<  tA4.coords[i]-tA3.coords[i] << std::endl;
            }
            
            return getTorsion3V(tR1, tR2, tR3);
        }
        else 
        {
            return 0.0;
        }
        
    }
    void Torsion::setValue()
    {
        
        if(atoms.size() == 4)
        {    
            itsValue = setValue(atoms[0], atoms[1], atoms[2], atoms[3]);
        }
        else
        {
            itsValue  = 0.0;
        }
    }
    
    
    REAL Torsion::getSigValue() const
    {
        return itsSigValue;
    }
    void Torsion::setSigValue(REAL tS)
    {
        itsSigValue = tS;
    }
    
    REAL Torsion::getValuePre() const
    {
        return itsValuePre;
    }
    void Torsion::setValuePre(REAL tP)
    {
        itsValuePre = tP;
    }
    
    REAL Torsion::getForceConst() const
    {
        return itsForceConst;
    }
    void Torsion::setForceConst(REAL tF)
    {
        itsForceConst = tF;
    }
    
    int Torsion::getPeriod() const
    {
        return itsPeriod;
    }
    void Torsion::setPeriod(int tP)
    {
        itsPeriod = tP;
    }
    
    // Another class for a torsion angle
    //Default constructor
    TorsionDict::TorsionDict():seriNum(-1),
            value(ZeroReal),
            sigValue(ZeroReal),
            valueST(ZeroReal),
            sigValueST(ZeroReal),
            period(1),
            id(NullString)
    {
    }
    
    // Copy constructor
    TorsionDict::TorsionDict(const TorsionDict& tTorsion):seriNum(tTorsion.seriNum),
            value(tTorsion.value),
            sigValue(tTorsion.sigValue),
            valueST(tTorsion.valueST),
            sigValueST(tTorsion.sigValueST),
            period(tTorsion.period),
            id(tTorsion.id)
    {
       
        for (std::vector<int>::const_iterator iAt = tTorsion.atoms.begin();
                iAt != tTorsion.atoms.end(); iAt++)
        {
            atoms.push_back(*iAt);
        }
        
        for (std::vector<AtomDict>::const_iterator iFA=tTorsion.fullAtoms.begin();
                iFA !=tTorsion.fullAtoms.end(); iFA++)
        {
            fullAtoms.push_back(*iFA);
        }
        
        for (std::vector<BondDict>::const_iterator iBo = tTorsion.bonds.begin();
                iBo != tTorsion.bonds.end(); iBo++)
        {
            bonds.push_back(*iBo);
        }
        
        for (std::vector<ID>::const_iterator iAt = tTorsion.atomCodClasses.begin();
                iAt != tTorsion.atomCodClasses.end(); iAt++)
        {
            atomCodClasses.push_back(*iAt);
        }
        
        for (std::vector<REAL>::const_iterator iAt = tTorsion.codTorsionValues.begin();
                iAt != tTorsion.codTorsionValues.end(); iAt++)
        {
            codTorsionValues.push_back(*iAt);
        }
    }
    
    // Constructor using atoms
    TorsionDict::TorsionDict(std::vector<AtomDict>& tAtoms):seriNum(-1),
            value(ZeroReal),
            sigValue(ZeroReal),
            valueST(ZeroReal),
            sigValueST(ZeroReal),
            period(1),
            id(NullString)
    {
        
        for (std::vector<AtomDict>::const_iterator iFA=tAtoms.begin();
                iFA !=tAtoms.end(); iFA++)
        {
            fullAtoms.push_back(*iFA);
        }
        
        for (std::vector<AtomDict>::iterator iAt = tAtoms.begin();
                iAt != tAtoms.end(); iAt++)
        {
            atomCodClasses.push_back(iAt->codClass);
        }
        
    }
    
    TorsionDict::~TorsionDict()
    {
    }
    
    extern int getTorsion(std::vector<TorsionDict> & tTors, 
                          int tAt1, int tAt2, int tAt3, int tAt4)
    {
        int tTor=-1;
        
        for (std::vector<TorsionDict>::iterator iT=tTors.begin();
                iT != tTors.end(); iT++)
        {
            //std::cout << "Current :" << std::endl
            //          << "atom1 " << iT->atoms[0] << "\t atom2 " << iT->atoms[1]
            //          << "\t atom 3 " << iT->atoms[2] << "\t atom 4 " 
            //          << iT->atoms[3] << std::endl;
           
            //if(std::find(iT->atoms.begin(), iT->atoms.end(), tAt1) !=iT->atoms.end()
            //   && std::find(iT->atoms.begin(), iT->atoms.end(), tAt2) !=iT->atoms.end()
            //   && std::find(iT->atoms.begin(), iT->atoms.end(), tAt3) !=iT->atoms.end()
            //   && std::find(iT->atoms.begin(), iT->atoms.end(), tAt4) !=iT->atoms.end())
            //{
            if ((iT->atoms[0]==tAt1 && iT->atoms[1]==tAt2 && iT->atoms[2]==tAt3 && iT->atoms[3]==tAt4)
                 || (iT->atoms[0]==tAt4 && iT->atoms[1]==tAt3 && iT->atoms[2]==tAt2 && iT->atoms[3]==tAt1))
            {
                tTor = iT->seriNum;
                // std::cout << "Torsion seriNumb " << iT->seriNum << std::endl;
                break;
            }
            
        }
        
        
        if (tTor ==-1)
        {
            std::cout << "can not find the target torsion "  << std::endl;
        }
          
        return tTor;
    }
    
 
    extern REAL getTorsion(std::vector<AtomDict> & tAtoms,
                           int iCur, int iNext, std::vector<int> tDoneSet)
    {
        
        int iPrev = tAtoms[iCur].tree["parent"][0];
        std::vector<int> tDoneCAs;
        for (std::vector<int>::iterator iC=tAtoms[iCur].connAtoms.begin();
                iC != tAtoms[iCur].connAtoms.end(); iC++)
        {
            if (*iC !=iNext && *iC !=iPrev 
                && std::find(tDoneSet.begin(), tDoneSet.end(), *iC) !=tDoneSet.end())
            {
                tDoneCAs.push_back(*iC);
            }
        }
        
        REAL retV =0.0;
        
        if ((int)tDoneCAs.size() ==0)
        {
            retV=tAtoms[iNext].treeTorsion;
        }
        
        return retV;
                
    }
    
    extern REAL getTorsion(AtomDict & tA1, 
                           AtomDict & tA2, 
                           AtomDict & tA3, 
                           AtomDict & tA4)
    {
        std::vector<REAL> tR1, tR2, tR3;
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
        if(tA1.coords.size() == 3 && tA2.coords.size() == 3 && 
                tA3.coords.size() ==3)
        {
            for (int i = 0; i < (int)tA1.coords.size(); i++)
            {
                tR1.push_back(tA2.coords[i]-tA1.coords[i]);
                //std::cout << "tR1 " << i+1 << "  "
                //        <<  tA2.coords[i]-tA1.coords[i] << std::endl;
                tR2.push_back(tA3.coords[i]-tA2.coords[i]);
                // std::cout << "tR2 " << i+1 << "  "
                //          <<  tA3.coords[i]-tA2.coords[i] << std::endl;
                tR3.push_back(tA4.coords[i]-tA3.coords[i]);
                // std::cout << "tR3 " << i+1 << "  "
                //          <<  tA4.coords[i]-tA3.coords[i] << std::endl;
            }
            
            return getTorsion3V(tR1, tR2, tR3);
        }
        else 
        {
            return 0.0;
        }
        
        return 0.0;
    }
    
}