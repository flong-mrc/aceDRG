
/* 
 * File:   restraintLists.cpp
 * Author: flong
 *
 * Created on September 22, 2011, 2:53 PM
 */

#include "restraintLists.h"

namespace LIBMOL 
{
    RestraintLists::RestraintLists()
    {      
    }
    RestraintLists::RestraintLists(RestraintLists& tRLists)
    {
        RestraintLists();
        
        if (tRLists.restrBondList.size()) 
        {
            for (std::vector<Bond>::iterator iT=tRLists.restrBondList.begin();
                    iT != tRLists.restrBondList.end(); iT++)
            {
                restrBondList.push_back((*iT));
            }
        }

        if (tRLists.restrAngleList.size())
        {          
            for (std::vector<Angle>::iterator iT=tRLists.restrAngleList.begin();
                    iT != tRLists.restrAngleList.end(); iT++)
            {
                restrAngleList.push_back((*iT));
            }
        }
        
        if (tRLists.restrTorsionList.size())
        {
            for (std::vector<Torsion>::iterator iT=tRLists.restrTorsionList.begin();
                    iT != tRLists.restrTorsionList.end(); iT++)
            {
                restrTorsionList.push_back((*iT));
            }
        }
        
        if (tRLists.restrChiralList.size())
        {
            for (std::vector<Chiral>::iterator iT= tRLists.restrChiralList.begin();
                    iT != tRLists.restrChiralList.end(); iT++)
            {
                restrChiralList.push_back((*iT));
            }
        }
        
        if (tRLists.restrPlaneList.size())
        {
            for (std::vector<Plane>::iterator iT=tRLists.restrPlaneList.begin();
                    iT != tRLists.restrPlaneList.end(); iT++)
            {
                restrPlaneList.push_back((*iT));
            }
        }
            
    }
    
    RestraintLists::~RestraintLists()
    {
        
        if (!restrBondList.empty())
        {
            restrBondList.clear();
        }
        
        if (!restrAngleList.empty())
        {
            restrAngleList.clear();
        }   
        
        if (!restrTorsionList.empty())
        {
            restrTorsionList.clear();
        }
        
        if (!restrChiralList.empty())
        {
            restrChiralList.clear();
        }
        
        if (!restrPlaneList.empty())
        {
            restrPlaneList.clear();
        }
        
    }
    
    //RestraintLists & RestraintLists::operator=(const RestraintLists & tR)
    //{    
    //}
}