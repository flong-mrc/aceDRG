/* 
 * File:   inputParams.cpp
 * Author: flong
 *
 * Created on November 1, 2012, 4:47 PM
 */

#include "inputParams.h"

namespace LIBMOL
{
    InputParams::InputParams() : aCifName(NullPoint),
            aMolSdfName(NullPoint),
            CodBondFileName(NullPoint),
            CodBondAngleFileName(NullPoint),
            userInputDir(NullPoint),
            userOutName(NullPoint),
            clibMonDir(NullString),
            monoRootName(NullString),
            workMode(1)
    {
    }
    
    InputParams::InputParams(int comArgc, char** comArgv): aCifName(NullPoint),
            aMolSdfName(NullPoint),
            CodBondFileName(NullPoint),
            CodBondAngleFileName(NullPoint),
            userInputDir(NullPoint),
            userOutName(NullPoint),
            clibMonDir(NullString),
            monoRootName(NullString),
            workMode(1)
    {
        std::string clibMonDir(std::getenv("CLIBD_MON"));
        if (clibMonDir.empty())
        {
            std::cerr << "You need to setup CCP4 suite first " << std::endl;
            exit(1);
        }
        
        if ( comArgc < 2)
        {
            printManu(1);
            exit(1);
        }
        
        parseInput(comArgc, comArgv);
        
    }
    
    InputParams::~InputParams()
    {
    }
    
    void InputParams::printManu(int nErr)
    {
        
    }
    
    void InputParams::parseInput(int comArgc, char** comArgv)
    {
        // do not need the following declarations, put here in case......
        //extern char *optarg;
        //extern int  optind, opterr, optopt;
        
        opterr = 0;
        int c, index;
        while ((c = getopt(comArgc, comArgv, "a:c:d:j:o:r:s:")) != -1)
        {
            switch (c)
            {
                case 'a':
                    CodBondAngleFileName = optarg;
                    break;
                case 'c':
                    aCifName =  optarg;
                    std::cout << "Input cif file is : " << aCifName << std::endl; 
                    workMode =1;
                    break;
                case 'd':
                    CodBondFileName = optarg;
                    break;
                case 'j':
                    userInputDir =  optarg;
                    std::cout << "The output root dir are in " 
                              <<  userInputDir << std::endl;
                    break;
                case 'o':
                    userOutName = optarg;
                    std::cout << "Output res file should be : " 
                              << userOutName << std::endl;
                    break;
                case 'r':
                    monoRootName = optarg;
                    std::cout << "Monomer root name is : " 
                        <<  monoRootName << std::endl;
                    break;
                case 's':
                    aMolSdfName  = optarg;
                    std::cout << "Input MolSdf file is : " << optarg 
                              << std::endl;
                    workMode = 2;
                    break;
                case '?':
                    if (std::isprint(optopt))
                    {
                        std::cerr <<  "command line arguments are not match at "
                                  << char(optopt) << std::endl;
                    }
                    else 
                    {
                        std::cerr << "Unknown option character: " << optopt 
                                  << std::endl;
                    }
                default:
                    continue;
            }
        }
        
        // Non-option arguments remained
        for (index = optind; index < comArgc; index++)
        {
            std::cout << "Non-option argument " << comArgv[index] 
                      << std::endl;
        }
    }
}