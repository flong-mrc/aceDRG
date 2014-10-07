 /* 
 * File:   atom.cpp
 * Author: flong
 *
 * last updated on October  1, 2011, 2:17 PM
 */


#include "atom.h"
#include "ring.h"

namespace LIBMOL
{
    // New Atom
    Atom::Atom():
            itsName(NullString),
            itsSeriNum(ZeroInt),
            itsElementType(NullString),
            itsID(NullString),
            itsResName(NullString),
            itsSeqNum(ZeroInt),
            itsInsCode("."),
            itsSegNum(ZeroInt),
            itsSegID(NullString),
            itsChainID(NullString),
            itsMass(ZeroReal),
            itsCharge(ZeroReal),
            itsPCharge(ZeroReal),
            itsRadius(ZeroReal),
            itsIonRadius(ZeroReal),
            itsAltLoc(NullString),
            itsOccup(ZeroReal),
            itsTempFact(ZeroReal),
            itsSigOcc(ZeroReal),
            itsSigTemp(ZeroReal)             
    {
        
        // allRestrLists = new RestraintLists();
    }

     Atom::Atom(const Atom & tAtom) 
     {
        // Atom();
         
         itsName        = tAtom.getName();
         itsSeriNum     = tAtom.getSeriNum();
         itsElementType = tAtom.getElementType();
         itsID          = tAtom.getID();
         itsResName     = tAtom.getResName();
         itsSeqNum      = tAtom.getSeqNum();
         itsChainID     = tAtom.getChainID();
         itsMass        = tAtom.getMass();
         itsInsCode     = tAtom.getInsCode();
         itsAltLoc      = tAtom.getAltLoc();
         
         for (int i = 0; i < (int)tAtom.Uxx.size(); i++)
         {
             Uxx.push_back(tAtom.Uxx[i]); 
         } 
         for (int i = 0; i <(int)tAtom.sigUxx.size(); i++)
         {
             sigUxx.push_back(tAtom.sigUxx[i]); 
         }         
         for (int i=0; i < (int)tAtom.coords.size(); i++)
         {
             coords.push_back(tAtom.coords[i]);
         }
         
         for (std::vector<REAL>::const_iterator iAtl=tAtom.altLoc.begin();
                 iAtl !=tAtom.altLoc.end(); iAtl++)
         {
             altLoc.push_back((*iAtl));
         }
     }
     
     Atom::~Atom()
     {
      /*   if(allRestrLists !=NULL)
         {
             delete allRestrLists;
             allRestrLists = NULL;
         }
       */
       
     }
    
    std::string Atom::getName() const
    {
        return itsName;
    }
    void Atom::setName(std::string tName)
    {
        itsName = tName;
    }
    
    int  Atom::getSeriNum() const
    {
        return itsSeriNum;
    }
    void Atom::setSeriNum(int tN)
    {
        itsSeriNum = tN; 
    }
    
    Element Atom::getElementType() const 
    {
        return itsElementType;
    }
    void Atom::setElementType(Element tElem)
    {
        itsElementType=tElem;
    }
        
    IDCode Atom::getID() const
    {
        return itsID;
    }
    void Atom::setIDCode(IDCode tID)
    {
        itsID=tID;
    }
        
    ResName Atom::getResName() const  
    {
        return itsResName;
    }
    void Atom::setResName(ResName tName)
    {
        itsResName=tName;
    }
        
    SeriNumber Atom::getSeqNum() const
    {
        return itsSeqNum;
    }
    void Atom::setSeqNum(SeriNumber tN)
    {
        itsSeqNum=tN;
    }
        
    InsCode Atom::getInsCode()  const
    {
        return itsInsCode;
    }
    void Atom::setInsCode(InsCode tID)
    {
        itsInsCode=tID;
    }
        
    SeriNumber Atom::getSegNum() const
    {
        return itsSegNum;
    }
    void  Atom::setSegNum(SeriNumber tN)
    {
        itsSegNum = tN;
    }
    
    IDCode Atom::getSegID()  const
    {
        return itsSegID;
    }
    void Atom::setSegID(IDCode tID)
    {
        itsSegID = tID;
    }
        
    ChainID Atom::getChainID() const
    {
        return itsChainID;
    }
    void Atom::setChainID(ChainID tID)
    {
        itsChainID = tID;
    }
      
    SeriNumber Atom::getModSeriNum() const  
    {
        return itsModNum;
    }
    void Atom::setModSeriNum(SeriNumber tN)
    {
        itsModNum = tN;
    }
        
    REAL Atom::getMass() const 
    {
        return itsMass;
    }
    void Atom::setMass(REAL tMass)
    {
        itsMass = tMass;
    }
        
    REAL Atom::getCharge()  const
    {
        return itsCharge;
    }
    void Atom::setCharge(REAL tC)
    {
        itsCharge = tC;
    }
        
    REAL Atom::getPartialCharge() const
    {
        return itsPCharge;
    }
    void Atom::setPartialCharge(REAL tC) 
    {
        itsPCharge = tC;
    }
        
    AltLoc Atom::getAltLoc() const
    {
        return itsAltLoc;
    }
    
    void Atom::setAltLoc(AltLoc tA)
    {
        itsAltLoc = tA;
    }
        
    REAL Atom::getOccup() const
    {
        return itsOccup;
    }
    void Atom::setOccup(REAL tO)
    {
        itsOccup = tO;
    }
        
    REAL Atom::getTempFact() const
    {
        return itsTempFact;
    }
    void Atom::setTempFact(REAL tT) 
    {
        itsTempFact = tT;
    }
        
    REAL Atom::getRadius() const
    {
        return itsRadius;
    }
    void Atom::setRadius(REAL tR)
    {
        itsRadius = tR;
    }
        
    REAL Atom::getIonRadius() const
    {
        return itsIonRadius;
    }
    void Atom::setIonRadius (REAL tIR)
    {
        itsIonRadius = tIR;
    }    
    
    // another kind of atom class 
    
    AtomDict::AtomDict() :seriNum(ZeroInt),
            resName(NullString),
            id(NullString),
            chemType(NullString),
            enerType(NullString),
            charge(ZeroReal),
            parCharge(ZeroReal),
            formalCharge(ZeroReal),
            radius(VDWCONST),
            ionRadius(VDWCONST),
            inChiralIdx(ZeroInt),
            chiralIdx(ZeroInt),
            bondingIdx(ZeroInt),
            isMetal(false),
            metalGeo(NullString),
            isoB(ZeroReal),
            ocp(1.0),
            symmMult(1),
            ccp4Type(NullString),
            cChemType(NullString),
            codClass(NullString),
            codMolIdx(ZeroInt),
            codCifName(NullString),
            codNBSymb(NullString),
            codNB2Symb(NullString),
            hashingValue(ZeroInt),
            coordExist(false),
            hasFreePi(false),
            isCChemTypeSet(false),
            isInPreCell(false),
            sId(NullString),
            symmOp(NullString),
            treeBond(ZeroReal),
            treeAngle(ZeroReal),
            treeTorsion(ZeroReal)
    {
        existProps["chemType"]   = -1;
        existProps["enerType"]   = -1;
        existProps["parCharge"]  = -1;
        existProps["x"]          = -1;
        existProps["y"]          = -1;
        existProps["z"]          = -1;
        int tDim=3;
        for (int i=0; i <  tDim; i++)
        {
            coords.push_back(0.0);
            fracCoords.push_back(0.0);
            forces.push_back(0.0);
        }
    }
    
    AtomDict::AtomDict(const AtomDict& tAtom) : seriNum(tAtom.seriNum), 
            resName(tAtom.resName),
            id(tAtom.id),
            chemType(tAtom.chemType),
            enerType(tAtom.enerType),
            charge(tAtom.charge),
            parCharge(tAtom.parCharge),
            formalCharge(tAtom.formalCharge),
            radius(tAtom.radius),
            ionRadius(tAtom.ionRadius),
            inChiralIdx(tAtom.inChiralIdx),
            chiralIdx(tAtom.chiralIdx),
            bondingIdx(tAtom.bondingIdx),
            isMetal(tAtom.isMetal),
            metalGeo(tAtom.metalGeo),
            isoB(tAtom.isoB),
            ocp(tAtom.ocp),
            symmMult(tAtom.symmMult),
            ccp4Type(tAtom.ccp4Type),
            cChemType(tAtom.cChemType),
            codClass(tAtom.codClass),
            codMolIdx(tAtom.codMolIdx),
            codCifName(tAtom.codCifName),
            codNBSymb(tAtom.codNBSymb),
            codNB2Symb(tAtom.codNB2Symb),
            hashingValue(tAtom.hashingValue),
            coordExist(tAtom.coordExist),
            hasFreePi(tAtom.hasFreePi),
            isCChemTypeSet(tAtom.isCChemTypeSet),
            isInPreCell(tAtom.isInPreCell),
            chiralChecked(tAtom.chiralChecked),
            sId(tAtom.sId),
            symmOp(tAtom.symmOp),
            treeBond(tAtom.treeBond),
            treeAngle(tAtom.treeAngle),
            treeTorsion(tAtom.treeTorsion)
    {
        for (std::map<std::string, int>::const_iterator iMa =tAtom.existProps.begin();
                iMa != tAtom.existProps.end(); iMa++)
        {
            existProps[iMa->first] = iMa->second;
        }
        
        for(std::vector<REAL>::const_iterator tX=tAtom.coords.begin();
                tX != tAtom.coords.end(); tX++)
        {
            coords.push_back((*tX));
        }
        for (int i=0; i < (int)tAtom.fracCoords.size(); i++)
         {
             fracCoords.push_back(tAtom.fracCoords[i]);
         }
        for(std::vector<REAL>::const_iterator tF=tAtom.forces.begin();
                tF != tAtom.forces.end(); tF++)
        {
            forces.push_back((*tF));
        }
        
        for (std::vector<int>::const_iterator iA=tAtom.connAtoms.begin();
                iA !=tAtom.connAtoms.end(); iA++)
        {
            connAtoms.push_back((*iA));
        }
        
        for (std::vector<int>::const_iterator iA=tAtom.connHAtoms.begin();
                iA !=tAtom.connHAtoms.end(); iA++)
        {
            connHAtoms.push_back((*iA));
        }
        
        for (std::vector<int>::const_iterator iA=tAtom.neighbAtoms.begin();
                iA != tAtom.neighbAtoms.end(); iA++)
        {
            neighbAtoms.push_back(*iA);
        }
        for (std::vector<int>::const_iterator iB = tAtom.inBonds.begin();
                iB != tAtom.inBonds.end(); iB++)
        {
            inBonds.push_back((*iB));
        }
        
        for (std::vector<AngleDict>::const_iterator iA=tAtom.inAngles.begin();
                iA != tAtom.inAngles.end(); iA++)
        {
            inAngles.push_back(*iA);
        }
        
        for (std::vector<int>::const_iterator iR = tAtom.inRings.begin();
                iR !=tAtom.inRings.end(); iR++)
        {
            inRings.push_back((*iR));
        }
        
        for(std::vector<int>::const_iterator iCh=tAtom.inChirals.begin();
                iCh!=tAtom.inChirals.end(); iCh++)
        {
            inChirals.push_back(*iCh);
        }
        for (std::vector<std::string>::const_iterator iNb = tAtom.nbRep.begin();
                iNb != tAtom.nbRep.end(); iNb++)
        {
            nbRep.push_back(*iNb);
        }
        for (std::map<std::string,int>::const_iterator iRiRep=tAtom.ringRep.begin();
                iRiRep != tAtom.ringRep.end(); iRiRep++)
        {
             ringRep.insert(std::pair<std::string,int>(iRiRep->first, iRiRep->second));
        }
        for (std::map<std::string,int>::const_iterator iRdRep=tAtom.ringRepBySeriNum.begin();
                iRdRep != tAtom.ringRepBySeriNum.end(); iRdRep++)
        {
             ringRepBySeriNum.insert(std::pair<std::string,int>(iRdRep->first, iRdRep->second));
        }
        for (std::vector<ID>::const_iterator iCod=tAtom.codClassV.begin();
                iCod !=tAtom.codClassV.end(); iCod++)
        {
            codClassV.push_back(*iCod);
        }
        for (std::map<std::string, std::vector<int> >::const_iterator iTr=tAtom.tree.begin();
                iTr != tAtom.tree.end(); iTr++)
        {
            for (std::vector<int>::const_iterator iAt = iTr->second.begin();
                    iAt != iTr->second.end(); iAt++)
            {
                tree[iTr->first].push_back(*iAt);
            }
        }
    }
    
    AtomDict::~AtomDict()
    {
    }
    
    int AtomDict::getNumAtomsWithin2stNB(std::vector<AtomDict> & tAllAtoms)
    {
        int tN=0;
        for (std::vector<int>::iterator iNB=connAtoms.begin();
                iNB !=connAtoms.end(); iNB++)
        {
            // allow repeat calculations in 2st NB
            tN +=(int)tAllAtoms[*iNB].connAtoms.size();    //should be (1+(int)tAllAtoms[*iNB].connAtoms-1)   
        }
        
        return tN;
    }
    
    /*
    void AtomDict::setCodClass()
    {
        codClass = "";
        
        codClass.append(id);
        outRingSec();
        
    }
    
    void AtomDict::outRingSec()
    {
        
        int numRings = (int)inRings.size();
        if (numRings)
        {
            for (int i = 0; i < numRings; i++)
            {
                int curRingSize = (int)inRings[i].atoms.size();
                std::string tSize = IntToStr(curRingSize);
                if (i==0)
                {
                    codClass.append("[");
                    if(numRings==1)
                    {
                        codClass.append(tSize + "]");
                    }
                    else
                    {
                        codClass.append(tSize + ",");
                    }
                }
                else
                {
                    if(i != numRings-1)
                    {
                        codClass.append(tSize + ",");
                    }
                    else
                    {
                        codClass.append(tSize + "]");
                    }
                } 
            }
        }
    }

     */
    
    int AtomDict::getMinRing() 
    {
        int rSize = 0;
        
        if ((int)ringRep.size() !=0)
        {
            int rMin = 7;
            for (std::map<ID, int>::iterator iR=ringRep.begin();
                    iR!=ringRep.end(); iR++)
            {
                if (iR->second < rMin)
                {
                    rMin = iR->second;
                }
            }
            rSize = rMin;
        }
        
        return rSize;
    }
    void AtomDict::setCodChemType()
    {
        // Change the second letter in element ID from upper to low case.
        if ((int)chemType.size()==1)
        {
            cChemType = chemType;
        }
        else if ((int)chemType.size() >1)
        {
            std::string tSubStr(chemType.substr(1));
            StrLower(tSubStr);
            cChemType = chemType.at(0)+ tSubStr;
        }
        
        // other properties for COD database
        
    }
    
    void AtomDict::fromCodClassToEnerType()
    {
        
    }
        
    void AtomDict::outNeighBAtoms()
    {
        
        // Check the atom
    }
 
   
    int AtomDict::atomPosition(std::vector<AtomDict>& tAtoms)
    {
        for (int i=0; i<(int)tAtoms.size(); i++)
        {
            if(tAtoms[i].id.compare(id)==0)
            {
                return i;
            }
        }
        
        return -1;
        
    }
   
}
