/* 
 * File:   codBondAndAngleGroups.h
 * Author: flong
 *
 * Created on August 30, 2012, 5:31 PM
 */

#ifndef CODBONDANDANGLEGROUPS_H
#define	CODBONDANDANGLEGROUPS_H
#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef FILE_H
#include "file.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef BOND_H
#include "bond.h"
#endif

#ifndef ANGLE_H
#include "angle.h"
#endif

namespace LIBMOL
{
    class AtomDict;
    class BondDict;
    class AngleDict;
    
    class CodBondAndAngleGroup
    {
    public:
        
        // default constructor
        CodBondAndAngleGroup();
        // the constructor using to file name
        CodBondAndAngleGroup(FileName               tBFName,
                             FileName               tAFName);
        
        // Destructor
        ~CodBondAndAngleGroup();
       
        void groupCodBonds(FileName       tBFName);
                          
        
        void groupCodAngles(FileName      tAFName);
        
        
        std::vector<BondDict>                    allDictBonds;
        std::map<int, std::map<int, std::map<ID, std::map<ID,  std::map<ID,
                 std::map<ID, std::map<ID, std::map<ID, int > > > > > > > >allDictBondsIdx;
        
        std::vector<AngleDict>       allDictAngles;
        std::map<int, std::map<int, std::map<int, 
        std::map<ID,  std::map<ID,  std::map<ID, 
        std::map<ID,  std::map<ID,  std::map<ID,
        std::map<ID,  std::map<ID,  std::map<ID,
        int > > > > > > > > > > > >  allDictAnglesIdx;
        
        
    };
}



#endif	/* CODBONDANDANGLEGROUPS_H */

