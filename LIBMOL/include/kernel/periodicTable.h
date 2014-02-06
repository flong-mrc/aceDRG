/* 
 * File:   periodcTable.h
 * Author: flong
 *
 * Created on July 29, 2012, 11:22 PM
 */

#ifndef PERIODICTABLE_H
#define	PERIODICTABLE_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

namespace LIBMOL
{
    
 
    class PeriodicTable
    {
    public :
        
        // Default constructor
        PeriodicTable();
        
        //Destructor
        ~PeriodicTable();
        
        std::map<ID,  std::map<ID, int> >                  elements;
        std::map<ID,  std::map<ID, REAL> >                 elemProps;
        std::map<int, std::string>                         matTypes;
        std::map<ID,  std::vector<std::pair<ID, int> > >   electronConf;
        std::map<ID, REAL>                                 electronegativities;
        std::map<ID, std::vector<int> >                    extraValences;
              
    };
    
}


#endif	/* PERIODCTABLE_H */

