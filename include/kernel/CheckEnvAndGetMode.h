/* 
 * File:   CheckEnvAndGetMode.h
 * Author: flong
 *
 * Created on January 10, 2013, 10:25 AM
 */

#ifndef CHECKENVANDGETMODE_H
#define	CHECKENVANDGETMODE_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef _GETOPT_H
#include <getopt.h>
#endif

#ifndef UTILITY_H
#include "utility.h"
#endif

#define noArgument 0
#define requiredArgument 1 
#define optionalArgument 2

namespace LIBMOL {

    class CheckEnvAndGetMode {
    public:
        
        // a global variable for acedrg tables 
        static const std::string aceDrgTabs;
        
        // Default constructor 
        CheckEnvAndGetMode();

        // Constructor take command-line arguments 
        CheckEnvAndGetMode(int numArg, char** argVars);

        // Default destructor
        ~CheckEnvAndGetMode();


        void parseArgs(int numArg, char** argVars);
        void setBandASiga();
        void SetWorkMode();
        void printManual();
        void printVarAndMode();
        
        int workMode;
        
        double                                   upperBondSig;
        double                                   lowBondSig;
        double                                   upperAngleSig;
        double                                   lowAngleSig;
        std::map<ID, ID> envVars;
        std::map<ID, ID> IOEntries;

    };

}


#endif	/* CHECKENVANDGETMODE_H */

