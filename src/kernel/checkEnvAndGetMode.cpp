/* 
 * File:   CheckEnvAndGetMode.cpp
 * Author: flong
 *
 * Created on January 10, 2013, 5:12 PM
 */

#include "CheckEnvAndGetMode.h"

namespace LIBMOL
{
    
    CheckEnvAndGetMode::CheckEnvAndGetMode():workMode(0)
    {
        // Currently CCP4 suite is the only requited thing 
        char * pClibdMon = std::getenv("CLIBD_MON");
        if (pClibdMon !=NULL)
        {
            if (envVars.find("CLIBD_MON") !=envVars.end())
            {
                envVars["CLIBD_MON"].clear();
            }
            envVars["CLIBD_MON"].append(pClibdMon);
        }
        else
        {
            std::cerr << "You need to setup CCP4 suite first " << std::endl;
            exit(1);
        }
        
        char * pLibMon = std::getenv("LIBMOL_ROOT");
        if (pLibMon !=NULL)
        {
           
            if (envVars.find("LIBMOL_ROOT") !=envVars.end())
            {
                envVars["LIBMOL_ROOT"].clear();
            }
            envVars["LIBMOL_ROOT"].append(pLibMon);
        }
        else
        {
            
            std::cerr << "You need to setup env variable $LIBMOL_ROOT first " << std::endl;
            exit(1);
        }
         
    }
    
    CheckEnvAndGetMode::CheckEnvAndGetMode(int numArg, char** ArgVars) :workMode(0)
    {
        // Currently CCP4 suite is the only requited thing 
        char * pClibdMon = std::getenv("CLIBD_MON");
        if (pClibdMon !=NULL)
        {
            if (envVars.find("CLIBD_MON") !=envVars.end())
            {
                envVars["CLIBD_MON"].clear();
            }
            envVars["CLIBD_MON"].append(pClibdMon) ;
        }
        else
        {
            std::cerr << "You need to setup CCP4 suite first " << std::endl;
            exit(1);
        }
       
        char * pLibMon = std::getenv("LIBMOL_ROOT");
        if (pLibMon !=NULL)
        {
            if (envVars.find("LIBMOL_ROOT") !=envVars.end())
            {
                envVars["LIBMOL_ROOT"].clear();
            }
            envVars["LIBMOL_ROOT"].append(pLibMon);
        }
        else
        {
            std::cerr << "You need to setup env variable $LIBMOL_ROOT first " << std::endl;
            exit(1);
        }
        
        // do not need the following declarations, put here in case......
        //extern char *optarg;
        // extern int optind, opterr, optopt;
        
        // opterr = 0;
           
        if (numArg  < 5)
        {
            // temporarily interface 
            std::cerr << "Command line syntax: " << std::endl
                      << "libmol -c CIFFILE -r monomer_name -o outputRestraintFile " 
                      << std::endl;
            exit(1);     
        }
        
        int c, index; 
        while ((c = getopt (numArg, ArgVars, "a:b:c:d:i:j:m:n:o:p:r:s:t:y:z:A:O:")) != -1)
        {
            switch (c)
            {
                case 'a':
                    IOEntries["CodBondAngleFileName"] = optarg;
                    break;
                case 'b':
                    IOEntries["inCifNameB"] = optarg;
                    break;    
                case 'c':
                    IOEntries["inCifName"]  =  optarg;
                    //std::cout << "Input cif file is : " << IOEntries["inCifName"]
                    //          << std::endl; 
                    break;     
                case 'd':
                    IOEntries["CodBondFileName"] = optarg;
                    break;
                case 'i':
                    IOEntries["inSmiName"] =  optarg;
                    // std::cout << "The input files are in " 
                    //          <<  IOEntries["inFileDir"] << std::endl;
                    break;    
                case 'j':
                    IOEntries["inFileDir"] =  optarg;
                    // std::cout << "The input files are in " 
                    //          <<  IOEntries["inFileDir"] << std::endl;
                    break;
                case 'm':
                    IOEntries["molGen"] =  optarg;
                    // std::cout << "The input files are in " 
                    //          <<  IOEntries["inFileDir"] << std::endl;
                    break;
                case 'n':
                    IOEntries["NBDepth"] =  optarg;
                    // std::cout << "The input files are in " 
                    //          <<  IOEntries["inFileDir"] << std::endl;
                    break;
                case 'o':
                    IOEntries["userOutName"] = optarg;
                    //std::cout << "Output res file should be : " 
                    //          << IOEntries["userOutName"] << std::endl;
                    break;
                case 'p':
                    IOEntries["inPdbName"] = optarg;
                    //std::cout << "Input PDB file is : " 
                    //          << IOEntries["inPdbName"] << std::endl;
                    break;
                case 'r':
                    IOEntries["monoRootName"] = optarg;
                    
                    //std::cout << "Monomer root name is : " 
                    //          << IOEntries["monoRootName"] << std::endl;
                    break;
                 case 's':
                    IOEntries["inSdfName"] = optarg;
                    //std::cout << "Monomer root name is : " 
                    //          << IOEntries["monoRootName"] << std::endl;
                    break; 
                 case 't':
                    IOEntries["tabGen"] = optarg;
                    std::cout << "Table generator mode :  " 
                              << IOEntries["tabGen"] << std::endl;
                    break;
                case 'y':
                    IOEntries["transCoords"] = optarg;
                    StrLower(IOEntries["transCoords"]);
                    break;                          
                case 'z':
                    IOEntries["transAngles"] = optarg;
                    StrLower(IOEntries["transAngles"]);
                    break;                    
                case 'A':
                    IOEntries["AtomTypeOutName"] = optarg;
                    //std::cout << "COD atom types are output to : " 
                    //          << IOEntries["AtomTypeOutName"] << std::endl;
                    break;
                case 'O':
                    IOEntries["NoGeoOpt"] = optarg;
                    StrLower(IOEntries["NoGeoOpt"]);
                    //std::cout << "COD atom types are output to : " 
                    //          << IOEntries["AtomTypeOutName"] << std::endl;
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
        for (index = optind; index < numArg; index++)
        {
            std::cout << "Non-option argument " << ArgVars[index] 
                      << std::endl;
        }
        
        if(IOEntries.find("NBDepth") == IOEntries.end())
        {
            IOEntries["NBDepth"] = "1";
        }
        
        std::map<ID,ID>::iterator iKeyFinder;
        iKeyFinder = IOEntries.find("monoRootName");
        if (iKeyFinder ==IOEntries.end() )
        {
            IOEntries["monoRootName"] ="UNL";
        }
        else
        {
            if (IOEntries["monoRootName"].empty())
            {
                IOEntries["monoRootName"] ="UNL";
            }
        }
        
        SetWorkMode();
        
    }
    
    CheckEnvAndGetMode::~CheckEnvAndGetMode()
    {
    }
    
    void CheckEnvAndGetMode::parseArgs(int numArg, char** argVars)
    {
        static int verbose_flag;
        
        static struct option longOptions[] =
        {
            /* These options set a flag. */
            {"verbose", noArgument,       &verbose_flag, 1},
            {"brief",   noArgument,       &verbose_flag, 0},
            /* These options don't set a flag.
               We distinguish them by their indices. 
            {"add",     noArgument,       0, 'a'}, */
            {"help",  noArgument,       0, 'h'}, 
            
            /* input files, keys are self-explained */
            {"cifin",  optionalArgument, 0, 'c'},
            {"pdbin",  optionalArgument, 0, 'p'},
            {"sdfin",  optionalArgument, 0, 'd'},
            {"smiin",  optionalArgument, 0, 'i'},  // the input smile string 
            /* output files, output pdb has the same name 
             * root as that of output cif  */
            {"cifout", optionalArgument, 0, 'o'},
            /* monomer name in output files */
            {"monout",  optionalArgument, 0, 'r'},
            {"codatm",  optionalArgument, 0, 'A'},
            {0, 0, 0, 0}
        };
        
        int c=0; 
        while (c !=-1)
        {
            int optIndex=0;
            c = getopt_long(numArg, argVars, "c:d:h:i:o:p:r:A:", longOptions, &optIndex);
            
            switch (c)
            {
                case 'h':
                    printManual();
                    break;
                case 'c':
                    IOEntries["inCifName"] =  optarg;
                    break;
                case 'p':
                    IOEntries["inPdbName"] =  optarg;
                    break;
                case 'd':
                    IOEntries["inSdfName"] =  optarg;
                    break; 
                case 'i':
                    IOEntries["inSmiStr"]  =  optarg;
                    break; 
                case 'o':
                    IOEntries["userOutName"] = optarg;
                    IOEntries["outCifName"]  =  optarg; // should be that
                    break; 
                case 'r':
                    IOEntries["monoRootName"] = optarg;
                    break; 
                case 'A':
                    IOEntries["AtomTypeOutName"] = optarg;
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
                    std::abort();
                    
            }
            
            
        }
    }
    
    void CheckEnvAndGetMode::printManual()
    {
        
    }
    void CheckEnvAndGetMode::SetWorkMode()
    {
        
        
        if ( IOEntries.find("inCifName")!=IOEntries.end())
        {
            if (!IOEntries["inCifName"].empty())
            {
                if (IOEntries.find("molGen") != IOEntries.end())
                {
                    workMode =32;
                }
                else if ( IOEntries.find("inPdbName")!=IOEntries.end())
                {
                    if (IOEntries.find("transCoords") != IOEntries.end())
                    {
                        if (!IOEntries["inPdbName"].empty() && IOEntries["transCoords"] =="y" )
                        {
                            // This mode read the input cif file and PDB file 
                            // replace the coordinates in the cif with the corresponding ones in the PDB
                            // which have been refined 
                            workMode = 21;
                            std::cout << "Coords transforming mode " << std::endl;
                        }
                        else
                        {
                            std::cout << "You need to enter a output PDB name "
                                     << std::endl;
                        }
                    }
                    if (IOEntries.find("transAngles") != IOEntries.end())
                    {
                        if (!IOEntries["inPdbName"].empty() && IOEntries["transAngles"] =="y" )
                        {
                            // This mode read the input cif file and PDB file 
                            // replace the angles in the cif with the ones calculated from PDB
                            workMode = 22;
                        }
                        else
                        {
                            std::cout << "You need to enter a output PDB name "
                                     << std::endl;
                        }
                    }
                }
                else
                {
                    if ( IOEntries.find("AtomTypeOutName")!=IOEntries.end())
                    {
                        if (!IOEntries["AtomTypeOutName"].empty())
                        {
                            // This mode outputs COD atom type only
                            workMode =41;
                        }
                        else
                        {
                            std::cout << "You need to enter the name of the output file containing COD atom type"
                                << std::endl;
                        }
                    }
                    else
                    {
                        // This is default restraint generation mode
                        workMode = 11;
                         
                        if (IOEntries.find("userOutName") ==IOEntries.end())
                        {
                            IOEntries["userOutName"] = "libmol.out"; 
                        }
                    }
                }
            }
        }
        else if (IOEntries.find("inCifNameB") !=IOEntries.end())
        {
            if (IOEntries.find("molGen") != IOEntries.end())
            {
                // This mode takes a general cif file and generate a whole molecule 
                // and unique bonds and angles within the molecule
                workMode = 31;
            }
            else if (IOEntries.find("tabGen") != IOEntries.end())
            {
                workMode = 33;
            }
        }
        else if (IOEntries.find("inSdfName") !=IOEntries.end())
        {
            if (IOEntries.find("monoRootName") != IOEntries.end())
            {
                workMode =13;
            }
            else if(IOEntries.find("AtomTypeOutName") !=IOEntries.end())
            {
                if (!IOEntries["AtomTypeOutName"].empty())
                {
                    // This mode outputs COD atom type only
                    workMode =43;
                }
            }
        }
        else if (IOEntries.find("inSmiName") !=IOEntries.end())
        {
            if (IOEntries.find("monoRootName") != IOEntries.end())
            {
                workMode =12;
            }
            else if(IOEntries.find("AtomTypeOutName") !=IOEntries.end())
            {
                if (!IOEntries["AtomTypeOutName"].empty())
                {
                    // This mode outputs COD atom type only
                    workMode =42;
                }
            }
        }
        else
        {
            std::cout << "what is the input file ? " << std::endl;
        }
        
        
        printVarAndMode();
        
    }
    
    void CheckEnvAndGetMode::printVarAndMode()
    {
        std::cout << "||=====================Job description====================================||"
                << std::endl;
        if (workMode==11 || workMode==12)
        {
            std::cout << "You will get an output cif file for a list of dictionary restraints"
                    << std::endl << "based on the input cif file" << std::endl
                    << "Your input cif file : " << IOEntries["inCifName"] << std::endl
                    << "Your monomer name : "   << IOEntries["monoRootName"]  << std::endl
                    << "The output dictionary file(cif) : " << IOEntries["userOutName"] 
                    << std::endl;
        }
        else if (workMode==12)
        {
            std::cout << "You will get an output cif file for a list of dictionary restraints"
                    << std::endl << "based on the input SMILE file" << std::endl
                    << "Your input SMILE file : " << IOEntries["inSmiStr"]      << std::endl
                    << "Your monomer name : "     << IOEntries["monoRootName"]  << std::endl
                    << "The output dictionary file(cif) : " << IOEntries["userOutName"] 
                    << std::endl;
        }
        else if (workMode==13)
        {
            std::cout << "You will get an output cif file for a list of dictionary restraints"
                    << std::endl << "based on the input MOL/SDF file"          << std::endl
                    << "Your input MOL/SDF file : " <<  IOEntries["inSdfName"] << std::endl
                    << "Your monomer name : "   << IOEntries["monoRootName"]   << std::endl
                    << "The output dictionary file(cif) : " << IOEntries["userOutName"] 
                    << std::endl;
        }
        else if (workMode==2)
        {
            std::cout << "Your job is get a full library description based on input cif and PDB files "
                    << std::endl
                    << "Your input cif file : " << IOEntries["inCifName"] << std::endl
                    << "Your input PDB file : " << IOEntries["inPdbName"] << std::endl
                    << "Your monomer name : "   << IOEntries["monoRootName"]  << std::endl
                    << "The output  dictionary file(cif) : " << IOEntries["userOutName"]
                    << std::endl;
        }
        else if (workMode==21)
        {
            std::cout << "The current job is coordinate transform in the dictionary file " 
                      << std::endl;
        }
        else if (workMode==31 || workMode==32)
        {
            std::cout << "Your job is to deduce values of bond and bond-angle"
                      << " and build molecules from your input cif file "
                      << std::endl
                      << "Your input cif file : " << IOEntries["inCifNameB"] << std::endl
                      << "Your monomer name : "   << IOEntries["monoRootName"]  << std::endl
                      << "The output dictionary file(cif) : " << IOEntries["userOutName"]
                      << std::endl;
        }   
        else if (workMode==41)
        {
            std::cout << "Your job is to look at atom types of COD format"
                    << std::endl << "for atoms in the input cif file" << std::endl
                    << "Your input cif file : " << IOEntries["inCifName"] << std::endl
                    << "Your monomer name : "   << IOEntries["monoRootName"]  << std::endl
                    << "The output atom type file(txt) : " << IOEntries["AtomTypeOutName"]
                    << std::endl;
        }
        else if (workMode==42)
        {
            std::cout << "Your job is to look at atom types of COD format"
                    << std::endl << "for atoms in the input mol/sdf file" << std::endl
                    << "Your input mol/sdf file : " << IOEntries["inSdfName"] << std::endl
                    << "Your monomer name : "   << IOEntries["monoRootName"]  << std::endl
                    << "The output atom type file(txt) : " << IOEntries["AtomTypeOutName"]
                    << std::endl;
        }
                  
        std::cout << "||========================================================================||"
                  << std::endl;
    }
}