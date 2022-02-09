/* 
 * File:   inputParams.h
 * Author: flong
 *
 * Created on November 1, 2012, 4:47 PM
 */

#ifndef INPUTPARAMS_H
#define	INPUTPARAMS_H

#ifndef _MSC_VER
#ifndef __UNISTD_
#include <unistd.h>
#endif
#endif

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef UTILITY_H
#include "utility.h"
#endif


namespace LIBMOL
{
    class InputParams 
    {
    public:
        // Default constructor
        InputParams();
        
        // Constructor using command line arguments 
        InputParams(int comArgc, char** comArgv);
        
        // Constructor using a xml file 
        InputParams(FileName aFName);
        
        // Destructor 
        ~InputParams();
        
        void printManu(int nErr);
        
        void parseInput(int comArgc, char** comArgv);
        
        
        FileName            aCifName;
        FileName            aMolSdfName;
        FileName            CodBondFileName;  
        FileName            CodBondAngleFileName;
        FileName            userInputDir;
        
        FileName            userOutName;
        
        std::string         clibMonDir;
        
        std::string         monoRootName;
        
        int                 workMode;  
    };
}


#endif	/* INPUTPARAMS_H */

