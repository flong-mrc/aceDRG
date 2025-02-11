/* 
 * 
 * File:   file.h
 * Author: flong
 *
 * Created on August 11, 2011, 4:06 PM
 */

#ifndef FILE_H
#define	FILE_H

#include <sys/types.h>
#include <sys/stat.h>

#ifndef KERNEL_H
#include "kernel.h"
#endif 

namespace LIBMOL
{
        enum fileType 
        {
            UNKNOWN = 0,
            isPATH,
            isRegularFile,
            isSymbolicLink,
            isSocketFile
        };
        
        struct fileState
        {
            bool   isOpen;
            bool   isReadable;
            bool   isWritable;
            bool   isExecutable;
        };
        
    class File : public  std::fstream 
    {
    public:
        File();
        
        
        // default constructor with file name and open mode
        File(std::string tName, std::ios::openmode tOpenMode);
        
        // destructor 
        virtual ~File()
        {
        }
        
        // some functions besides those in std::fstream
        

        // getName just get the base name without path
        Name                    getName() throw();    
        void                    setName(std::string) throw();
        std::string             getFullPath() throw();
        void                    setFullPath() throw();
        Name                    getFileExtName() throw(); 
        Size                    getSize();
        std::ios::openmode      getOpenMode();
        fileType                getType();
         
        
        bool                    copyFile(std::string tTarget,
                                         std::size_t tBuff=4096) throw();
        bool                    moveFile(std::string tTarget) throw();
        bool                    removeFile() throw();
        
    private :
        
        Name                    itsName;
        std::string             itsFullPath;
        Name                    itsExt;
        Size                    itsSize;
        std::ios::openmode      itsOpenMode;
        fileType                itsType;
        fileState               itsState;
         
    };
}


#endif	/* FILE_H */

