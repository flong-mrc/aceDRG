/* 
 * File:   ExtraRestrDictFile.h
 * Author: flong
 *
 * Created on September 13, 2011, 5:02 PM
 */

#ifndef EXTRARESTRDICTFILE_H
#define	EXTRARESTRDICTFILE_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef UTILITY_H
#include "utility.h"
#endif

#ifndef FILE_H
#include "file.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef ATOMASSEMBLY_H
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

#ifndef RESTRAINTLISTS_H
#include "restraintLists.h"
#endif


namespace LIBMOL
{
    class Atom;
    class Residue;
    class ModRes;
    class Chain;
    class Model;

    class Bond;
    class Angle;
    class Torsion;
    class Chiral;
    class Plane;
    
    class Helix;
    class Sheet;
    class Strand;
    class Turn;
    
    class SSBond;
    class Link;
    
    class RestraintLists;
    
    enum RestrInfoType
    {
        // currently, only two of them
        DICT_ADDTIONAL,
        DICT_PAIRS,
        BOND,
        ANGLE,
        TORSION,
        CHIRAL,
        PLANE      
    };
    
    struct RestrOneMonoInfo
    {
        std::string monomerName;
        std::string groupName;
    };
    // better inherited from class 'RestraitLists' (do it later)
    class RestraitListSet
    {  
    public:
        // 
        RestraitListSet();
        RestraitListSet(const RestraitListSet & tP);
        ~RestraitListSet();
        
        std::string                   label;   
        std::string                   groupLabel;
        std::vector<RestrOneMonoInfo> monoSet;
        std::vector<Bond>             restrBondList; 
        std::vector<Angle>            restrAngleList;
        std::vector<Torsion>          restrTorsionList;
        std::vector<Chiral>           restrChiralList;
        std::vector<Plane>            restrPlaneList;
    private:
        
    };
    
    class ExtraRestrDictFile : public File
    {
    public :
        
        // default constructor
        ExtraRestrDictFile();
        
        // copy constructor 
        ExtraRestrDictFile(ExtraRestrDictFile & tF);
        ExtraRestrDictFile(Name tF);
        ExtraRestrDictFile(FileName tF);
        // destructor 
        ~ExtraRestrDictFile();
        
        void setupSystem();
        
        RestrInfoType getRecordType(ALine tLine);
        bool setAddRestrDefInfo(ALine  tLine);
        bool setPairRestrDefInfo(ALine  tLine);
        
        void setInfoType();
        bool setOneBondRestr(ALine           tLine);
        bool setOneAngleRestr(ALine          tLine);
        bool setOneTorsionRestr(ALine        tLine);
        bool setOneChiralRestr(ALine         tLine);
        bool setOnePlaneRestr(ALine          tLine);
        
        
        RestrInfoType      getRestrInfoType();
        void               setRestrInfoType(RestrInfoType tR);
             
        
        std::ifstream      inFile;

        std::vector<RestraitListSet>     allRestrLists;
        
    private :
        
        // name, open mode etc and their associated access methods
        // are defined in 'File' class 
        Name                    itsName;
        std::string             itsFullPath;
        Name                    itsExt;
        Size                    itsSize;
        std::ios::openmode      itsOpenMode;
        fileType                itsType;
        fileState               itsState;
       
        RestrInfoType           itsInfoType;
        RestraitListSet *       itsCurRestrListSet;
        RestrOneMonoInfo  *     itsCurOneMono;
        
    };
    
}

#endif	/* FREERESTRFILE_H */

