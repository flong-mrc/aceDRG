/* 
 * File:   MonomerLib.h
 * Author: flong
 *
 * Created on September 27, 2011, 12:53 PM
 */

#ifndef MONOMERLIB_H
#define	MONOMERLIB_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef FILE_H
#include "file.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef LIBG_ATOMASSEMBLY_H
#include "atomAssembly.h"
#endif

#ifndef RESIDUE_H
#include "residue.h"
#endif

#ifndef CHAIN_H
#include "chain.h"
#endif

#ifndef LIBG_MODEL_H
#include "libgmodel.h"
#endif

#ifndef SSBOND_H
#include "ssbond.h"
#endif

#ifndef LIBG_LINK_H
#include "libglink.h"
#endif

#ifndef CRYSTINFO_H
#include "crystInfo.h"
#endif

#ifndef SECONDARYSTRUCTURES_H
#include "secondaryStructures.h"
#endif

namespace LIBMOL
{
    class MonomerLib 
    {
    public :
        
        // default constructor 
        MonomerLib();
        
        // MonomerLib constructor via a dictionary file name
        MonomerLib(FileName tF);
        
        // destructor
        ~MonomerLib();
        
        // exception classes 
        class setMonomerException
        {
        };
        
        FileName     getMonoListName();
        void         setMonoListName(FileName tF);
        
        void         setupSystem(FileName  tF);
        
        // public member for quick access 
        std::map<std::string, std::string>    monomerGroups;
        
        
    private :
        
        FileName               itsMonoListName;
        
    };
}

#endif	/* MONOMERLIB_H */

