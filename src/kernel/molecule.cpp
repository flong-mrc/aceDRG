
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
            hasCoords(true),
            validated(true)
    {
    }
    
    Molecule::Molecule(const Molecule & tMol) : seriNum(tMol.seriNum),
            id(tMol.id),
            formula(tMol.formula),
            hasCoords(tMol.hasCoords),
            validated(tMol.validated)
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
}