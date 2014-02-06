
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
            hasCoords(true)
    {
    }
    
    Molecule::Molecule(const Molecule & tMol) : seriNum(tMol.seriNum),
            id(tMol.id),
            hasCoords(tMol.hasCoords)
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
}