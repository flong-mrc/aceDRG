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
            altId(NullString),
            chemType(NullString),
            enerType(NullString),
            charge(ZeroReal),
            parCharge(ZeroReal),
            formalCharge(ZeroReal),
            formalChargeI(ZeroInt),
            radius(VDWCONST),
            ionRadius(VDWCONST),
            inChiralIdx(ZeroInt),
            chiralIdx(ZeroInt),
            bondingIdx(ZeroInt),
            hybrid(NullString),
            matType(0),
            isMetal(false),
            metalGeo(NullString),
            isoB(ZeroReal),
            ocp(1.0),
            symmMult(1),
            excessElec(0),
            ccp4Type(NullString),
            cChemType(NullString),
            codClass(NullString),
            codMolIdx(ZeroInt),
            codCifName(NullString),
            codAtmRoot(NullString),
            codNBSymb(NullString),
            codNB2Symb(NullString),
            codNB3Symb(NullString),
            codAtmMain(NullString),
            codNB1NB2_SP(NullString),
            codNB1NB2_ExElec(NullString),
            hashingValue(ZeroInt),
            coordExist(false),
            numPi(ZeroInt),
            isCChemTypeSet(false),
            isInPreCell(false),
            chiralChecked(false),
            isInAromRing(false),
            isInSP2Ring(false),
            fromCalc(false),
            sId(NullString),
            symmOp(NullString),
            fromOrig(-1),
            sMolType(NullString),
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

        baseRingProp["size"]  ="";
        baseRingProp["aroma"] ="";

        formalCharge = charge;
    }

    AtomDict::AtomDict(const AtomDict& tAtom) : seriNum(tAtom.seriNum),
            resName(tAtom.resName),
            id(tAtom.id),
            altId(tAtom.altId),
            chemType(tAtom.chemType),
            enerType(tAtom.enerType),
            charge(tAtom.charge),
            parCharge(tAtom.parCharge),
            formalCharge(tAtom.formalCharge),
            formalChargeI(tAtom.formalChargeI),
            radius(tAtom.radius),
            ionRadius(tAtom.ionRadius),
            inChiralIdx(tAtom.inChiralIdx),
            chiralIdx(tAtom.chiralIdx),
            bondingIdx(tAtom.bondingIdx),
            hybrid(tAtom.hybrid),
            matType(tAtom.matType),
            isMetal(tAtom.isMetal),
            metalGeo(tAtom.metalGeo),
            isoB(tAtom.isoB),
            ocp(tAtom.ocp),
            symmMult(tAtom.symmMult),
            excessElec(tAtom.excessElec),
            ccp4Type(tAtom.ccp4Type),
            cChemType(tAtom.cChemType),
            codClass(tAtom.codClass),
            codMolIdx(tAtom.codMolIdx),
            codCifName(tAtom.codCifName),
            codAtmRoot(tAtom.codAtmRoot),
            codNBSymb(tAtom.codNBSymb),
            codNB2Symb(tAtom.codNB2Symb),
            codNB3Symb(tAtom.codNB3Symb),
            codAtmMain(tAtom.codAtmMain),
            codNB1NB2_SP(tAtom.codNB1NB2_SP),
            codNB1NB2_ExElec(tAtom.codNB1NB2_ExElec),
            hashingValue(tAtom.hashingValue),
            coordExist(tAtom.coordExist),
            numPi(tAtom.numPi),
            isCChemTypeSet(tAtom.isCChemTypeSet),
            isInPreCell(tAtom.isInPreCell),
            chiralChecked(tAtom.chiralChecked),
            isInAromRing(tAtom.isInAromRing),
            isInSP2Ring(tAtom.isInSP2Ring),
            fromCalc(tAtom.fromCalc),
            sId(tAtom.sId),
            symmOp(tAtom.symmOp),
            fromOrig(tAtom.fromOrig),
            sMolType(tAtom.sMolType),
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

        for (std::vector<int>::const_iterator iA=tAtom.connMAtoms.begin();
                iA !=tAtom.connMAtoms.end(); iA++)
        {
            connMAtoms.push_back((*iA));
        }

        for (std::vector<int>::const_iterator iA=tAtom.neighbAtoms.begin();
                iA != tAtom.neighbAtoms.end(); iA++)
        {
            neighbAtoms.push_back(*iA);
        }

        for (std::map<int, REAL>::const_iterator iM=tAtom.NBAtomMap.begin();
                iM != tAtom.NBAtomMap.end(); iM++)
        {
            NBAtomMap[iM->first]=iM->second;
        }

        for (std::vector<int>::const_iterator iB = tAtom.inBonds.begin();
                iB != tAtom.inBonds.end(); iB++)
        {
            inBonds.push_back((*iB));
        }

        for (std::vector<REAL>::const_iterator iBL=tAtom.bondLengths.begin();
                iBL != tAtom.bondLengths.end(); iBL++)
        {
            bondLengths.push_back(*iBL);
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

        for (std::vector<RingDict>::const_iterator
                iRF=tAtom.inRingsFull.begin();
                iRF != tAtom.inRingsFull.end(); iRF++)
        {
            inRingsFull.push_back(*iRF);
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

        for (std::map<std::string,std::string>::const_iterator iRiRep=tAtom.ringRepS.begin();
                iRiRep != tAtom.ringRepS.end(); iRiRep++)
        {
             ringRepS.insert(std::pair<std::string,std::string>(iRiRep->first, iRiRep->second));
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

        for (std::map<std::string, std::string>::const_iterator iMRP=tAtom.baseRingProp.begin();
                iMRP != tAtom.baseRingProp.end(); iMRP++)
        {
            baseRingProp[iMRP->first] = iMRP->second;
        }

        for (std::vector<std::string>::const_iterator
             iFMType=tAtom.formType.begin();
             iFMType !=tAtom.formType.end(); iFMType++)
        {
            formType.push_back(*iFMType);
        }

        formalCharge = charge;
    }

    AtomDict::AtomDict(const Atom & tAtom):seriNum(tAtom.getSeriNum()-1),
                                resName(tAtom.getResName()),
                                id(tAtom.getName()),
                                altId(NullString),
                                chemType(tAtom.getElementType()),
                                enerType(NullString),
                                charge(tAtom.getPartialCharge()),
                                parCharge(ZeroReal),
                                formalCharge(ZeroReal),
                                formalChargeI(ZeroInt),
                                radius(VDWCONST),
                                ionRadius(VDWCONST),
                                inChiralIdx(ZeroInt),
                                chiralIdx(ZeroInt),
                                bondingIdx(ZeroInt),
                                hybrid(NullString),
                                matType(0),
                                isMetal(false),
                                metalGeo(NullString),
                                isoB(ZeroReal),
                                ocp(1.0),
                                symmMult(1),
                                excessElec(0),
                                ccp4Type(NullString),
                                cChemType(NullString),
                                codClass(NullString),
                                codMolIdx(ZeroInt),
                                codCifName(NullString),
                                codAtmRoot(NullString),
                                codNBSymb(NullString),
                                codNB2Symb(NullString),
                                codNB3Symb(NullString),
                                codAtmMain(NullString),
                                codNB1NB2_SP(NullString),
                                codNB1NB2_ExElec(NullString),
                                hashingValue(ZeroInt),
                                coordExist(false),
                                numPi(ZeroInt),
                                isCChemTypeSet(false),
                                isInPreCell(false),
                                chiralChecked(false),
                                isInAromRing(false),
                                isInSP2Ring(false),
                                fromCalc(false),
                                sId(NullString),
                                symmOp(NullString),
                                fromOrig(-1),
                                sMolType(NullString),
                                treeBond(ZeroReal),
                                treeAngle(ZeroReal),
                                treeTorsion(ZeroReal)
    {
        for(std::vector<REAL>::const_iterator tX=tAtom.coords.begin();
                tX != tAtom.coords.end(); tX++)
        {
            coords.push_back((*tX));
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

    int AtomDict::getNum1stNbHave2edNb(std::vector<AtomDict>& tAtoms)
    {
        int t12=0;
        for (std::vector<int>::iterator iNB=connAtoms.begin();
                iNB != connAtoms.end(); iNB++)
        {
            if (tAtoms[*iNB].connAtoms.size() >= 2)
            {
                t12++;
            }
        }

        return t12;

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

    int AtomDict::getMinRing2()
    {
        int rSize = 0;

        if (codClass.size() !=0)
        {
            std::vector<ID> tmpStrs1, tmpStrs2, tmpStrs3, tmpStrs4;
            StrTokenize(codClass, tmpStrs1, '(');
            if (tmpStrs1.size() > 0)
            {
                if (tmpStrs1[0].find("[") !=std::string::npos)
                {
                    StrTokenize(tmpStrs1[0],  tmpStrs2, '[');
                    if (tmpStrs2.size() > 1)
                    {
                        if (tmpStrs2[1].find(",") !=std::string::npos)
                        {
                            StrTokenize(tmpStrs2[1],  tmpStrs3, ',');
                            if (tmpStrs3.size() >0)
                            {
                                if (tmpStrs3[0].find("x") !=std::string::npos)
                                {
                                    StrTokenize(tmpStrs3[0],  tmpStrs4, 'x');
                                    if (tmpStrs4.size() >1)
                                    {
                                        rSize = StrToInt(tmpStrs4[1]);
                                    }
                                }
                                else
                                {
                                    rSize = StrToInt(tmpStrs3[0]);
                                }
                            }
                        }
                        else
                        {
                            StrTokenize(tmpStrs2[1],  tmpStrs3, ']');
                            if (tmpStrs3.size() >0)
                            {
                                if (tmpStrs3[0].find('x') !=std::string::npos)
                                {
                                    StrTokenize(tmpStrs3[0],  tmpStrs4, 'x');
                                    if (tmpStrs4.size() >1)
                                    {
                                        rSize = StrToInt(tmpStrs4[1]);
                                    }
                                }
                                else
                                {
                                    rSize = StrToInt(tmpStrs3[0]);
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            rSize = getMinRing();
        }

        return rSize;
    }

    void AtomDict::setBaseRingProps()
    {
        if (codClass.size() !=0)
        {
            ID tClass = TrimSpaces(codClass);
            std::vector<ID> tSecs;
            StrTokenize(tClass, tSecs, '(');
            if (tSecs.size() > 0)
            {
                if (tSecs[0].find("[") !=std::string::npos)
                {
                    std::vector<ID> tRS;
                    StrTokenize(tSecs[0], tRS, '[');
                    if (tRS.size() > 1)
                    {
                        // std::cout << "atom type " << codClass << std::endl;
                        // std::cout << "tRS[1] " << tRS[1]  << std::endl;
                        if (tRS[1].find("a") !=std::string::npos)
                        {
                            baseRingProp["aroma"] = "y";
                        }
                        // std::cout << "ring prop " << baseRingProp["aroma"] << std::endl;
                        cleanChar(tRS[1], ']');

                        if(tRS[1].find(",") !=std::string::npos)
                        {
                            std::vector<ID> tMinR;
                            StrTokenize(tRS[1], tMinR, ',');
                            if (tMinR.size() >0)
                            {
                                if (tMinR[0].find("x")==std::string::npos)
                                {
                                    if (tMinR[0].find("a")==std::string::npos)
                                    {
                                        baseRingProp["size"] = tMinR[0];
                                    }
                                    else
                                    {
                                        std::vector<ID> tOS;
                                        StrTokenize(tMinR[0], tOS, 'a');
                                        baseRingProp["size"] = tOS[0];
                                    }
                                }
                            }
                        }
                        else
                        {
                            if (tRS[1].find("x")==std::string::npos)
                            {
                                if (tRS[1].find("a")==std::string::npos)
                                {
                                    baseRingProp["size"] = tRS[1];
                                }
                                else
                                {
                                    std::vector<ID> tOS;
                                    StrTokenize(tRS[1], tOS, 'a');
                                    baseRingProp["size"] = tOS[0];
                                }
                            }
                        }
                    }
                }
            }
        }
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


    metalCluster::metalCluster():metSeril(-1),
                                 formu(NullString),
                                 coordGeoID(NullString)
    {
    }

    metalCluster::metalCluster(const metalCluster & tMC):metSeril(tMC.metSeril),
                                                        formu(tMC.formu),
                                                        coordGeoID(tMC.coordGeoID)
    {
        for (std::vector<int>::const_iterator iLA=tMC.ligandSerilSet.begin();
                iLA !=tMC.ligandSerilSet.end(); iLA++)
        {
            ligandSerilSet.push_back(*iLA);
        }

        for (std::map<int, std::map<int, std::vector<int> > >::const_iterator
                  iLF=tMC.ligandNBs.begin();
                  iLF!=tMC.ligandNBs.end(); iLF++)
        {
            for (std::map<int, std::vector<int> >::const_iterator
                    iOneLF=iLF->second.begin();
                    iOneLF !=iLF->second.end(); iOneLF++)
            {
                for (std::vector<int>::const_iterator
                     iNB = iOneLF->second.begin();
                     iNB != iOneLF->second.end(); iNB++)
                {
                    ligandNBs[iLF->first][iOneLF->first].push_back(*iNB);
                }
            }
        }

        for (std::map<int, std::string >::const_iterator
                iFo=tMC.ligandForma.begin();
                iFo != tMC.ligandForma.end(); iFo++)
        {
            ligandForma[iFo->first] = iFo->second;
        }

        for (std::map<int, std::map<int, REAL> >::const_iterator
                iBM=tMC.uniqBondsMap.begin();
                iBM != tMC.uniqBondsMap.end(); iBM++)
        {
            for (std::map<int, REAL>::const_iterator iBV=iBM->second.begin();
                    iBV !=iBM->second.end(); iBV++)
            {
                uniqBondsMap[iBM->first][iBV->first] = iBV->second;
            }
        }

        for (std::map<int, std::map<int, std::map<int, REAL> > >::const_iterator
                iAngM=tMC.uniqAngsMap.begin(); iAngM !=tMC.uniqAngsMap.end();
                iAngM++)
        {
            for (std::map<int, std::map<int, REAL> >::const_iterator
                    iAng1=iAngM->second.begin();
                    iAng1 != iAngM->second.end(); iAng1++)
            {
                for (std::map<int, REAL>::const_iterator
                      iAng2=iAng1->second.begin();
                      iAng2 !=iAng1->second.end(); iAng2++)
                {
                    uniqAngsMap[iAngM->first][iAng1->first][iAng2->first]
                           = iAng2->second;
                }
            }
        }

        for (std::map<int, std::map<int, std::map<int, REAL> > >::const_iterator
            iAngM=tMC.secondNBAngsMap.begin(); iAngM !=tMC.secondNBAngsMap.end();
            iAngM++)
        {
            for (std::map<int, std::map<int, REAL> >::const_iterator
                    iAng1=iAngM->second.begin();
                    iAng1 != iAngM->second.end(); iAng1++)
            {
                for (std::map<int, REAL>::const_iterator
                      iAng2=iAng1->second.begin();
                      iAng2 !=iAng1->second.end(); iAng2++)
                {
                    secondNBAngsMap[iAngM->first][iAng1->first][iAng2->first]
                           = iAng2->second;
                }
            }
        }

        allRings.clear();
        for (std::vector<RingDict>::const_iterator iR=tMC.allRings.begin();
              iR != tMC.allRings.end(); iR++)
        {
            allRings.push_back(*iR);
        }

    }

    metalCluster::~metalCluster()
    {

    }

    void metalCluster::setMetClusterFormu(std::vector<AtomDict>& tAtoms)
    {
        // The metal cluster formula.
        formu = setOneAtomFormu(tAtoms, metSeril);
        std::cout << "Metal cluster formular " << formu << std::endl;

        // Formula for each ligand atoms
        for (std::vector<int>::iterator iL=ligandSerilSet.begin();
                iL !=ligandSerilSet.end(); iL++)
        {
            ligandForma[*iL]= setOneAtomFormu(tAtoms, *iL);
            std::cout << "NB for ligand atom " << tAtoms[*iL].id
                      << " is " << ligandForma[*iL] << std::endl;
        }
    }


    std::string metalCluster::setOneAtomFormu(std::vector<AtomDict>& tAtoms,
                                                 int tCenSerial)
    {

        std::string aReturn = "";
        int tOrig;
        if (!tAtoms[tCenSerial].isInPreCell)
        {
            tOrig = tAtoms[tCenSerial].fromOrig;
        }
        else
        {
            tOrig = tAtoms[tCenSerial].seriNum;
        }

        aReturn += (tAtoms[tOrig].chemType + ":");

        std::map<std::string, std::vector<int> > tmpLF;

        for (std::vector<int>::iterator iLS=tAtoms[tOrig].connAtoms.begin();
                iLS !=tAtoms[tOrig].connAtoms.end(); iLS++)
        {
            tmpLF[tAtoms[*iLS].chemType].push_back(*iLS);
        }

        std::string aLast = "";
        std::list<std::string> tContainer;
        for (std::map<std::string, std::vector<int> >::iterator iT=tmpLF.begin();
                iT != tmpLF.end(); iT++)
        {
            if (iT->first.compare("H")!=0)
            {
                tContainer.push_back("{" + iT->first
                                     + IntToStr(iT->second.size())
                                     + "}");
            }
            else
            {
                // H atoms should not exist in the cluster, put in the last
                // for debug
                aLast = "{" + iT->first + IntToStr(iT->second.size())
                        + "}";
            }
        }

        tContainer.sort(compareNoCase2);

        for (std::list<std::string>::iterator iComp=tContainer.begin();
                iComp != tContainer.end(); iComp++)
        {
             aReturn+=*iComp;
        }

        if(aLast.size() !=0)
        {
            aReturn +=aLast;
            if (tCenSerial==metSeril)
            {
                std::cout << "A bug: metal atom " << tAtoms[metSeril].id
                      << " of serial number " << metSeril
                      << " connects directly to H atoms " << aLast
                      << ". Please check!" << std::endl;
            }

        }

        // How to record the charge of the cluster ?

        return aReturn;

    }

    void metalCluster::buildBondAndAngleMap(std::vector<AtomDict> & tAtoms,
                                       std::vector<CrystInfo>::iterator tCryst)
    {
        std::cout << "Build the metal cluster centered at atom "
                  << tAtoms[metSeril].id << " of serial number "
                  << tAtoms[metSeril].seriNum << std::endl;
        std::vector<REAL>  bVList, aVList, a2VList;

        for (std::vector<int>::iterator iNB = ligandSerilSet.begin();
                    iNB != ligandSerilSet.end(); iNB++)
        {


            REAL rD = getBondLenFromFracCoords(tAtoms[metSeril].fracCoords,
                                               tAtoms[*iNB].fracCoords,
                                tCryst->itsCell->a, tCryst->itsCell->b,
                                tCryst->itsCell->c, tCryst->itsCell->alpha,
                                tCryst->itsCell->beta, tCryst->itsCell->gamma);

            if (!inVectABS(bVList,  rD, 0.00001))
            {
                std::cout << "neighbor atom " << tAtoms[*iNB].id
                      << " of serial number " << tAtoms[*iNB].seriNum << std::endl;
                std::cout << "Distance " << rD << std::endl;
                // filter the exact same bond values caused by the symm atoms
                uniqBondsMap[metSeril][*iNB] = rD;
                bVList.push_back(rD);
            }

            // build angle map for Metal--NB1--NB22
            for (std::vector<int>::iterator iNB22=ligandNBs[metSeril][*iNB].begin();
                 iNB22!=ligandNBs[metSeril][*iNB].end(); iNB22++)
            {
                REAL rA22 = getAngleValueFromFracCoords(tAtoms[*iNB],
                                        tAtoms[metSeril], tAtoms[*iNB22],
                                        tCryst->itsCell->a, tCryst->itsCell->b,
                                        tCryst->itsCell->c, tCryst->itsCell->alpha,
                                        tCryst->itsCell->beta, tCryst->itsCell->gamma);
                REAL tRA22 = rA22*PID180;
                if (!inVectAllABS(a2VList, tRA22, 0.0001))
                {
                    secondNBAngsMap[metSeril][*iNB][*iNB22] = tRA22;
                    a2VList.push_back(tRA22);
                }

            }

            for (std::vector<int>::iterator iNB2 = ligandSerilSet.begin();
                                iNB2 != ligandSerilSet.end(); iNB2++)
            {
                if (*iNB2 > *iNB)
                {
                    //std::cout   << " the other NB " << tAtoms[*iNB2].id
                    //            << " of serial " << tAtoms[*iNB2].seriNum
                    //            << std::endl;

                    REAL rA = getAngleValueFromFracCoords(tAtoms[metSeril],
                                        tAtoms[*iNB], tAtoms[*iNB2],
                                        tCryst->itsCell->a, tCryst->itsCell->b,
                                        tCryst->itsCell->c, tCryst->itsCell->alpha,
                                        tCryst->itsCell->beta, tCryst->itsCell->gamma);
                    // std::cout << "Angle " << rA * PID180 << std::endl;
                    REAL tRA = rA*PID180;
                    if (!inVectAllABS(aVList, tRA, 0.0001))
                    {
                        uniqAngsMap[metSeril][*iNB][*iNB2] = tRA;
                        aVList.push_back(tRA);
                    }
                }
            }
        }
        std::cout << "The whole cluster is done." << std::endl;
    }

    void metalCluster::projectCoordsToUnitSph(std::vector<AtomDict>& tAtoms)
    {
        unitCoords.clear();

        std::vector<REAL>   tZero;
        tZero.push_back(0.0);
        tZero.push_back(0.0);
        tZero.push_back(0.0);
        unitCoords.push_back(tZero);

        for (std::vector<int>::iterator iNB = ligandSerilSet.begin();
                    iNB != ligandSerilSet.end(); iNB++)
        {
            std::vector<REAL> aDiffVec;
            aDiffVec.clear();
            if (tAtoms[*iNB].coords.size()==tAtoms[metSeril].coords.size())
            for (unsigned i=0; i < tAtoms[*iNB].coords.size(); i++)
            {
                aDiffVec.push_back(tAtoms[*iNB].coords[i] -tAtoms[metSeril].coords[i]);
            }

            normalizeV(aDiffVec);

            unitCoords.push_back(aDiffVec);
        }
    }

    void metalCluster::setOrgRingsInMetalCluster(std::vector<AtomDict> & tAtoms,
                                                 int tCenSerial)
    {
    }

    extern int getAtom(std::string              tId,
                       int                      tSeri,
                       std::vector<AtomDict> &  tAtoms)
    {
        int lFind= -1;

        for (std::vector<AtomDict>::iterator iAt=tAtoms.begin();
                iAt !=tAtoms.end(); iAt++)
        {
            if (iAt->seriNum==tSeri && iAt->id.compare(tId)==0)
            {
                lFind = iAt->seriNum;
                break;
            }
        }

        return lFind;
    }


    extern void setAtomsAltId(std::vector<AtomDict> & tAtoms)
    {
        bool lDub = false;

        std::map<std::string, std::vector<std::string> >  atomNameDistri;

        for (std::vector<AtomDict>::iterator iAt=tAtoms.begin();
             iAt !=tAtoms.end(); iAt++)
        {
            atomNameDistri[iAt->id].push_back(iAt->id);
            if (atomNameDistri[iAt->id].size() > 1)
            {
                lDub= true;
                break;
            }
        }

        if (lDub)
        {

            std::map<std::string, int> elemDistr;

            for (std::vector<AtomDict>::iterator iAt=tAtoms.begin();
                 iAt !=tAtoms.end(); iAt++)
            {
                if (elemDistr.find(iAt->chemType)==elemDistr.end())
                {
                    elemDistr[iAt->chemType] =1;
                }
                else
                {
                    elemDistr[iAt->chemType] +=1;
                }
            }

            for (std::vector<AtomDict>::iterator iAt=tAtoms.begin();
                 iAt !=tAtoms.end(); iAt++)
            {
                iAt->altId = iAt->chemType + IntToStr(elemDistr[iAt->chemType]);
                elemDistr[iAt->chemType]--;
            }
        }
        else
        {
            for (std::vector<AtomDict>::iterator iAt=tAtoms.begin();
                 iAt !=tAtoms.end(); iAt++)
            {
                iAt->altId = iAt->id;
            }
        }
        // Check
        //std::cout << "Now atoms are renamed " << std::endl;

        //for (std::vector<AtomDict>::iterator iAt=tAtoms.begin();
        //        iAt !=tAtoms.end(); iAt++)
        //{
        //    std::cout << "Original name : " << iAt->id
        //              << "  New name : " << iAt->altId << std::endl;
        //}

    }

    extern void setMetalConnAtm(std::vector<AtomDict> & tAtoms)
    {

        for (std::vector<AtomDict>::iterator iAt = tAtoms.begin();
                iAt != tAtoms.end(); iAt++)
        {
            iAt->connMAtoms.clear();
            for (std::vector<int>::iterator iNB = iAt->connAtoms.begin();
                         iNB != iAt->connAtoms.end(); iNB++)
            {
                if (tAtoms[*iNB].isMetal)
                {
                    iAt->connMAtoms.push_back(*iNB);
                }
            }
        }

    }
}
