/* 
 * File:   PDBHeaderSec.h
 * Author: flong
 *
 * Created on September 9, 2011, 3:57 PM
 */

#ifndef PDBHEADERSEC_H
#define	PDBHEADERSEC_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

namespace LIBMOL
{

        struct HeaderInfo 
        {
          Name                recordName;
          std::string         classification;
          Date                depositDate;
          ID                  PDBId;
        };
        
        struct ObsoleteInfo 
        {
           Name                recordName;
           ID                  PDBId;
           ID                  replacePDBId;
           Date                replaceDate;
        }; 
        
        struct TitleInfo
        {
            Name               recordName;
            int                contIndicator;
            std::string        title; 
        };
        
        struct SplitInfo 
        {
            Name               recordName;
            ID      *          complexParts;                
        };
        
        struct CaveatInfo
        {
            Name               recordName;
            ID                 entryID;
            std::string        comment;
        };
        
        struct CompndInfo
        {
            Name               recordName;
            std::string   *    compound;
        };
        
        struct SourceInfo
        {
            Name               recordName;
            std::string   *    srcName;
            
        };
        
        struct KeywdsInfo 
        {
            Name               recordName;
            std::string   *    keyWords;
        };
 
                // structure that stores the header records in a PDB file
        struct HeadSec
        {
            struct HeaderInfo          *    HEADER;
            struct ObsoleteInfo        *    OBSLTE;
            struct TitleInfo           *    TITLE;
            struct SplitInfo           *    SPLIT;
            struct CaveatInfo          *    CAVEAT;
            struct CompndInfo          *    COMPND; 
            struct SourceInfo          *    SOURCE;
            struct KeywdsInfo          *    KEYWDS;
        };
        
}

#endif	/* PDBHEADERSEC_H */

