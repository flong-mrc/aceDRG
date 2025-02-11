/* 
 * File:   link.h
 * Author: flong
 *
 * Created on September 6, 2011, 11:44 PM
 */

#include "libglink.h"

namespace LIBMOL
{
    Link::Link():itsName(NullString),
            itsID(NullString),
            itsSeriNum(ZeroInt),
            itsLength(ZeroReal)
    {
        atoms.empty();
    }
    
    Link::Link(const Link & tLink): itsName(tLink.getName()),
            itsID(tLink.getID()),
            itsSeriNum(tLink.getSeriNum()),
            itsLength(ZeroReal)
    {
        for (std::vector<Atom>::const_iterator iA=tLink.atoms.begin(); 
                iA != tLink.atoms.end(); iA++)
        {
            atoms.push_back((*iA));
        }
    }
    
    Link::~Link()
    {
        if (atoms.size())
        {
            atoms.clear();
        }
    }
    
     Name Link::getName() const
    {
            return itsName;
    }
    void Link::setName(Name tNa)
    {
        itsName = tNa;
    }
        
    ID Link::getID() const
    {
        return itsID;
    }
    void Link::setID(ID tID)
    {
        itsID = tID;
    }
        
    SeriNumber Link::getSeriNum() const
    {
        return itsSeriNum;
    }
    void Link::setSeriNum(SeriNumber tSer)
    {
        itsSeriNum = tSer;
    }
       
  
    void  Link::setLength()
    {
    
        itsLength = ZeroReal;
        
        if (atoms.size() == 2)
        {
            if(atoms[0].coords.size() && atoms[1].coords.size())
            {
                for (int i =0; i <(int)atoms[0].coords.size(); i++)
                {
                    itsLength +=((atoms[0].coords[i]- atoms[1].coords[i])
                            *(atoms[0].coords[i]- atoms[1].coords[i]));
                }
            }
            
            itsLength = std::sqrt(itsLength);       
        }
    }        
    
}