/* 
 * File:   libgmodel.cpp
 * Author: flong
 *
 * Created on September 13, 2011, 2:22 PM
 */

#include "libgmodel.h"

namespace LIBMOL
{
    Model::Model()
    {
    }
     
    Model::Model(const Model& tM)
    {
        for (int i=0; i < (int)tM.getNumOfChains(); i++)
        {
            addOneChain(tM.chains[i]);
        }
        
        for (int i=0; i < (int)tM.getNumOfHelices(); i++)
        {
            addOneHelix(tM.helices[i]);
        }
        
        for (int i=0; i < (int)tM.getNumOfSheets(); i++)
        {
            addOneSheet(tM.sheets[i]);
        }        
    }
    
    Model::~Model()
    {
        if(!chains.empty())
        {
            chains.clear();
        }
        if(!links.empty())
        {
            links.clear();
        }
        if(!helices.empty())
        {
            helices.clear();
        }
        if(!sheets.empty())
        {
            sheets.clear();
        }
    }
    

        
    ID Model::getID() const
    {
        return itsID;
    }
    void Model::setID(ID tID)
    {
        itsID = tID;
    }
        
    SeriNumber Model::getSeriNum() const
    {
        return itsSeriNum;
    }
    void Model::setSeriNum(SeriNumber tSer)
    {
        itsSeriNum = tSer;
    }
          
    // for chains
    Size Model::getNumOfChains() const
    {
        return chains.size();
    }
    
    Chain & Model::getOneChain(SeriNumber tI) 
    {
        return chains[tI];
    }
    
    void Model::addOneChain(const Chain & tC)
    {
        chains.push_back(tC);
    }
    
    bool Model::deleteOneChain(SeriNumber tI)
    {
        if (tI < (int)chains.size())
        {
            chains.erase(chains.begin()+tI-1);
            return true;
        }
        
        return false;
    }
    // for Links
    
    Size Model::getNumOfLinks() const
    {
        return links.size();
    }
    
    Link & Model::getOneLink(SeriNumber tN)
    {
        return links[tN];
    }
    
    void Model::addOneLink(const Link& tL)
    {
        links.push_back(tL);
    }
    
    bool Model::deleteOneLink(SeriNumber tN)
    {
        if (tN < (int)links.size())
        {
            links.erase(links.begin()+tN-1);
            return true;
        }
        
        return false;
    }
    
    // for helices
    
    Size Model::getNumOfHelices() const
    {
        return helices.size();
    }
    
    Helix & Model::getOneHelix(SeriNumber tI) 
    {
        return helices[tI];  
    }
    
    void Model::addOneHelix(const Helix& tH)
    {
        helices.push_back(tH);
    }
    
    void Model::deleteOneHelix(SeriNumber tI)
    {
        helices.erase(helices.begin()+tI-1);
    }
    
    // For sheets
    
    Size Model::getNumOfSheets() const
    {
        return sheets.size();
    }
    
    Sheet & Model::getOneSheet(SeriNumber tI) 
    {
       return sheets[tI];
    }
    
    void Model::addOneSheet(const Sheet & tS)
    {
        sheets.push_back(tS);
    }
    
    bool Model::deleteOneSheet(SeriNumber tN)
    {
        if (tN < (int)sheets.size())
        {
            sheets.erase(sheets.begin()+tN-1);
            return true;
        }
        
        return false;  
    }

    
}