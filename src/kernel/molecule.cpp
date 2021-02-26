
/* 
 * File:   molecule.cpp
 * Author: flong
 *
 * Created on November 1, 2012, 1:49 PM
 */

#include "molecule.h"

namespace LIBMOL
{
    
    Molecule::Molecule() : seriNum(-1),
            id(NullString),
            formula(NullString), 
            sumExcessElecs(ZeroInt),
            sumCharges(ZeroInt),
            hasCoords(true),
            validated(true), 
            isInf(false),
            stateChanged(false)
    {
    }
    
    Molecule::Molecule(const Molecule & tMol) : seriNum(tMol.seriNum),
            id(tMol.id),
            formula(tMol.formula),
            hasCoords(tMol.hasCoords),
            sumExcessElecs(tMol.sumExcessElecs),
            sumCharges(tMol.sumCharges),
            validated(tMol.validated),
            isInf(tMol.isInf),
            stateChanged(tMol.stateChanged)
    {
        for(std::vector<AtomDict>::const_iterator iA=tMol.atoms.begin();
                iA !=tMol.atoms.end(); iA++)
        {
            atoms.push_back(*iA);        
        }
        
        for (std::vector<BondDict>::const_iterator iB=tMol.bonds.begin();
                iB != tMol.bonds.end(); iB++)
        {
            bonds.push_back(*iB);
        }
        
        for (std::vector<BondDict>::const_iterator iB=tMol.allBonds.begin();
                iB != tMol.allBonds.end(); iB++)
        {
            allBonds.push_back(*iB);
        }
        
        for (std::vector<AngleDict>::const_iterator iAn=tMol.angles.begin();
                iAn !=tMol.angles.end(); iAn++)
        {
            angles.push_back(*iAn);
        }
        
        for (std::vector<TorsionDict>::const_iterator iT=tMol.torsions.begin();
                iT !=tMol.torsions.end(); iT++)
        {
            torsions.push_back(*iT);
        }
        
        for (std::vector<RingDict>::const_iterator iR=tMol.rings.begin();
                iR !=tMol.rings.end(); iR++)
        {
            rings.push_back(*iR);
        }
        
        for (std::vector<ChiralDict>::const_iterator iC=tMol.chirals.begin();
                iC != tMol.chirals.end(); iC++)
        {
            chirals.push_back(*iC);
        }     
        
        for (std::vector<AtomDict>::const_iterator iHA=tMol.extraHAtoms.begin();
                iHA !=tMol.extraHAtoms.end(); iHA++)
        {
            extraHAtoms.push_back(*iHA);
        }
        for (std::map<ID, std::vector<std::string> >::const_iterator iProp=tMol.propData.begin();
                iProp != tMol.propData.end(); iProp++)
        {
            for (std::vector<std::string>::const_iterator iVal=iProp->second.begin();
                    iVal !=iProp->second.end(); iVal++)
            {
                propData[iProp->first].push_back(*iVal);
            }
        }
        
        for (std::vector<std::string>::const_iterator iComm=tMol.comments.begin();
                iComm != tMol.comments.end(); iComm++)
        {
            comments.push_back(*iComm);
        }
        
    }
    
    Molecule::~Molecule()
    {
        
    }
    
    
    void Molecule::setAtomCartCoordFromFracCoord(
                                      std::vector<CrystInfo>::iterator tCryst)
    {
        for (std::vector<AtomDict>::iterator iAtm=atoms.begin();
                iAtm != atoms.end(); iAtm++)
        {
            if (iAtm->fracCoords.size()==3)
            {
                FractToOrtho(iAtm->fracCoords, iAtm->coords, 
                             tCryst->itsCell->a, tCryst->itsCell->b,
                             tCryst->itsCell->c, tCryst->itsCell->alpha,
                             tCryst->itsCell->beta, tCryst->itsCell->gamma);
            }
        }
    }
    void Molecule::setFormula()
    {
        std::map<ID, int> fmap;
        
        for (std::vector<AtomDict>::iterator iAt=atoms.begin();
                iAt != atoms.end(); iAt++)
        {
            if (iAt->formalCharge==0.0)
            {
                if (fmap.find(iAt->chemType) == fmap.end())
                {
                    fmap[iAt->chemType] = 1;
                }
                else
                {
                    fmap[iAt->chemType]++;
                }
            }
            else
            {
                ID tReal, tInt, tKey;
                tReal = RealToStr(iAt->formalCharge);
                
                
                if (tReal.find('.') !=tReal.npos)
                {
                    std::vector<ID> tBuf;
                    StrTokenize(tReal, tBuf, '.');
                    if (tBuf.size()!=0)
                    {
                        if(tBuf[0] !="0")
                        {
                            tInt = tBuf[0];
                        }
                    }
                    
                    if (tInt.size() !=0)
                    {
                        tKey = iAt->chemType + "[" + tInt + "]";
                    }
                    else
                    {
                        tKey = iAt->chemType;
                    }    
                }
                else if (tReal.size() !=0)
                {
                    if (tReal.substr(0,1).find("0") == std::string::npos)
                    {
                        tKey = iAt->chemType + "[" + tReal+ "]";
                    }
                    else
                    {
                        tKey = iAt->chemType;
                    }  
                }
                else
                {
                    tKey = iAt->chemType;
                }
                
                if (fmap.find(tKey) == fmap.end())
                {
                    fmap[tKey] = 1;
                }
                else
                {
                    fmap[tKey]++;
                } 
            }
        }
        
        std::list<std::string> tAllIds;
                            
        for (std::map<ID, int>::iterator iM=fmap.begin();
                iM != fmap.end(); iM++)
        {
            ID aID = "{" + iM->first + IntToStr(iM->second) + "}";
            
            tAllIds.push_back(aID);
                               
        }
        
        tAllIds.sort(compareNoCase);
                            
        formula = "";
        for (std::list<std::string>::iterator iAI =tAllIds.begin();
                iAI != tAllIds.end(); iAI++)
        {
            formula.append(*iAI);
        }  
    }
    
    void Molecule::calcSumExcessElecs()
    {
        sumExcessElecs = sumExElectrons(atoms);
    }
    
    void Molecule::calcSumCharges()
    {
        sumCharges = 0;
        
        // Unify charges and formal-charges on atoms
        for (std::vector<AtomDict>::iterator iAt=atoms.begin();
                iAt != atoms.end(); iAt++)
        {
            if (iAt->formalCharge !=0.0)
            {
                iAt->charge = iAt->formalCharge;
            }
            else
            {
                if (iAt->charge !=0)
                {
                    iAt->formalCharge = iAt->charge;
                }
            }
        }
        
        for (std::vector<AtomDict>::iterator iAt=atoms.begin();
                iAt != atoms.end(); iAt++)
        {
            sumCharges +=(iAt->formalCharge);
        }
    }
    
    void Molecule::setAtomFormTypes()
    {
        for (std::vector<AtomDict>::iterator iAtm=atoms.begin();
                    iAtm !=atoms.end(); iAtm++)
        {
            if (iAtm->chemType.compare("H")!=0)
            {
                iAtm->formType.clear();
                iAtm->formType.push_back(iAtm->chemType);
                if (iAtm->isInAromRing)
                {
                    std::cout << "Atom " << iAtm->id 
                              << " is at least in one aromatic ring " 
                              << std::endl;
                }
                std::string sAll="", s1="", s2="", s3="";
                s1 = iAtm->chemType + "_sp" + IntToStr(iAtm->bondingIdx);
                if (iAtm->isInAromRing)
                {
                    s2 = "_arom";
                }
            
                
                //if (iAtm->isInDelocBond)
                // {
                //    s3 = "_deloc";
                //}
                
            
                if (s2.size()>0 && s3.size() > 0)
                {
                    sAll = s1 + s2 +s3;
                }
                else if (s2.size() > 0)
                {
                    sAll = s1 + s2;
                }
                else if (s3.size() >0)
                {
                    sAll = s1 + s3;
                }
                else
                {
                    sAll = s1;
                }
            
                iAtm->formType.push_back(sAll);
            }
        }
        
        for (std::vector<AtomDict>::iterator iAtm=atoms.begin();
                    iAtm !=atoms.end(); iAtm++)
        {
            if (iAtm->chemType.compare("H")==0)
            {
                iAtm->formType.clear();
                iAtm->formType.push_back(iAtm->chemType);
                
                std::string aS;
                
                if (iAtm->connAtoms.size() >0)
                {
                    if (atoms[iAtm->connAtoms[0]].formType.size()==2)
                    {
                        aS =iAtm->chemType + "_"
                            + atoms[iAtm->connAtoms[0]].formType[1];
                    }
                    else if (atoms[iAtm->connAtoms[0]].formType.size()==1)
                    {
                        aS =iAtm->chemType + "_"
                            + atoms[iAtm->connAtoms[0]].formType[0];
                    }
                    else
                    {
                        aS = iAtm->chemType + "_" 
                            +  atoms[iAtm->connAtoms[0]].chemType;
                    }
                }
                else
                {
                    aS = iAtm->chemType + "_X";
                }
                
                iAtm->formType.push_back(aS);
            }
        }
        

    }
}