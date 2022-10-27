/*
 * 
 * File:   CIFFile.cpp
 * Author: flong
 *
 * last updated on  Jan 24, 2012, 05:36 PM
*/

#include "DictCifFile.h"
#ifdef _MSC_VER
#include <ciso646>
#endif

namespace LIBMOL
{
    
    
    GenCifFile::GenCifFile() : curBlockLine(ZeroInt),
            lErr(false),
            hasCoords(false),
            hasH(false),
            hasMetal(false),
            symmOprOK(true),
            notPowder(true),
            notNeuD(true),
            resolOK(true),
            RFactorOK(false),
            colidOK(true),
            hasOcpLab(true),
            hasHeavyCalcAtoms(false),
            nonCheck(false),
            checkR(true),
            checkResol(true),
            RTHRESHOLD_U(0.10),
            RESOLTHRESHOLD_U(0.842),
            itsCurAtomSeriNum(ZeroInt),
            itsCurAtom(NullPoint),
            itsCurBlock(""),
            itsCurCryst(NullPoint)
    {
    }
    
    GenCifFile::GenCifFile(Name tFname, 
                           std::ios_base::openmode tOpenMode) :
                           curBlockLine(ZeroInt),
                           lErr(false),
                           hasCoords(false),
                           hasH(false),
                           hasMetal(false),
                           symmOprOK(true),
                           notPowder(true),
                           notNeuD(true),
                           resolOK(true),
                           RFactorOK(false),
                           colidOK(true),
                           hasOcpLab(true),
                           hasHeavyCalcAtoms(false),
                           nonCheck(false),
                           checkR(true),
                           checkResol(true),
                           RTHRESHOLD_U(0.10),
                           RESOLTHRESHOLD_U(0.842),
                           itsCurAtomSeriNum(ZeroInt),  
                           itsCurAtom(NullPoint),
                           itsCurBlock(""),
                           itsCurCryst(NullPoint)
    {
        
       
        if (tOpenMode == std::ios::in)
        {
            inFile.open(tFname.c_str(), tOpenMode);
            if (inFile.is_open())
            {
                itsCurAtom    = new AtomDict();
                setupSystem2();
            }
            else
            {
                std::cout << tFname << " Can not be opened for reading " << std::endl;
                errMsg.push_back("Cif file Can not be opened for reading ");
                lErr = true;
            }
        }
        else
        {
            outFile.open(tFname.c_str(), tOpenMode);
        }
        
    }
    
    GenCifFile::GenCifFile(Name tFname, 
                           std::string    tParaFName,
                           std::ios_base::openmode tOpenMode) :
                           curBlockLine(ZeroInt),
                           lErr(false),
                           hasCoords(false),
                           hasH(false),
                           hasMetal(false),
                           symmOprOK(true),
                           notPowder(true),
                           notNeuD(true),
                           resolOK(true),
                           RFactorOK(false),
                           colidOK(true),
                           hasOcpLab(true),
                           hasHeavyCalcAtoms(false),
                           nonCheck(false),
                           checkR(true),
                           checkResol(true),
                           RTHRESHOLD_U(0.10),
                           RESOLTHRESHOLD_U(0.842),
                           itsCurAtomSeriNum(ZeroInt),  
                           itsCurAtom(NullPoint),
                           itsCurBlock(""),
                           itsCurCryst(NullPoint)
    {
         
        std::ifstream  aFilterParaF(tParaFName.c_str());
        if (aFilterParaF.is_open())
        {
            
            std::vector<std::string>     allLines;
            std::string                  aRecord;
            while(!aFilterParaF.eof())
            {   
                std::getline(aFilterParaF, aRecord);
                aRecord = TrimSpaces(aRecord);
                allLines.push_back(aRecord);
            }
            aFilterParaF.close();
            if (allLines.size() > 0)
            {
                setAllCrits(allLines);
            }
        }
        
        /*
        else
        {
            std::string aErrMsg = tParaFName + " Can not be opened for reading";
            std::cout << aErrMsg << std::endl;
            errMsg.push_back(aErrMsg + "\n");
            lErr = true;
            
        }
         */
        if (tOpenMode == std::ios::in)
        {
            inFile.open(tFname.c_str(), tOpenMode);
            if (inFile.is_open())
            {
                itsCurAtom    = new AtomDict();
                setupSystem2();
            }
            else
            {
                std::cout << tFname << " Can not be opened for reading " << std::endl;
                errMsg.push_back("Cif file Can not be opened for reading ");
                lErr = true;
            }
        }
        else
        {
            outFile.open(tFname.c_str(), tOpenMode);
        }
        
    }
    
    GenCifFile::GenCifFile(FileName                tFname,
                           std::ios::openmode      tOpenMode=std::ios::in ) :
                           curBlockLine(ZeroInt),
                           lErr(false),
                           hasCoords(false),
                           hasH(false),
                           hasMetal(false),
                           symmOprOK(true),
                           notPowder(true),
                           notNeuD(true),
                           resolOK(true),
                           RFactorOK(false),
                           colidOK(true),
                           hasOcpLab(true),
                           hasHeavyCalcAtoms(false),
                           nonCheck(false),
                           checkR(true),
                           checkResol(true),
                           RTHRESHOLD_U(0.10),
                           RESOLTHRESHOLD_U(0.842),
                           itsCurAtomSeriNum(ZeroInt),
                           itsCurAtom(NullPoint),
                           itsCurCryst(NullPoint)
    {   
        
        if (tOpenMode == std::ios::in)
        {
            inFile.open(tFname, tOpenMode);
            if (inFile.is_open())
            {
                itsCurAtom    = new AtomDict();
            
                setupSystem();
            }
            else
            {
                std::cout << tFname << " Can not be opened for reading " << std::endl; 
            }
        }
        else
        {
            outFile.open(tFname, tOpenMode);
        }        
    }
    
    
    
    GenCifFile::~GenCifFile()
    {
        if(inFile.is_open())
        {
            inFile.close();
        }
        
        if(outFile.is_open())
        {
            outFile.close();
        }
        if (itsCurAtom)
        {
            delete itsCurAtom; 
            itsCurAtom = NULL;
        }
        if(itsCurCryst)
        {
            delete itsCurCryst;
            itsCurCryst = NULL;
        }
    }
    
    void GenCifFile::setupSystem()
    {   
        
        if (inFile.is_open() )
        { 
            std::vector<std::vector<std::string> > tBlocs;
            // make sure
            itsCurBlock = "";
            std::string tRecord="";
            
            std::vector<std::string>     tAllLines;
            std::vector<std::string>     tBlocLines;
            
            while(!inFile.eof())
            {   
                std::getline(inFile, tRecord);
                tRecord = TrimSpaces(tRecord);
                tAllLines.push_back(tRecord);
                // std::cout <<  tRecord << std::endl;
                std::vector<std::string> tBuf;
                std::vector<std::string> tBuf_t;
                // only use a few blocks in the cif file
                
                
                if (tRecord.find("loop_") !=std::string::npos 
                    || tRecord.find("data_") !=std::string::npos)
                {
                    if (!tBlocLines.empty())
                    {
                        tBlocs.push_back(tBlocLines);
                        tBlocLines.clear();
                    }
                    
                    itsCurBlock = "loop";
                }
                else if (itsCurBlock=="loop")
                {
                    tBlocLines.push_back(tRecord);
                }
            }
            
            if (!tBlocLines.empty())
            {
                tBlocs.push_back(tBlocLines);
                tBlocLines.clear();
            }
            inFile.close();
            
            
            std::cout << "Number of data blocks in the cif file is " << tBlocs.size() << std::endl;

            
            if ((int)tBlocs.size()!=0)
            {    
                for (std::vector<std::vector<std::string> >::iterator iBs=
                        tBlocs.begin(); iBs !=tBlocs.end(); iBs++)
                {
                    
                    if((int)iBs->size() !=0)
                    {
                        for (std::vector<std::string>::iterator iBl=iBs->begin();
                               iBl != iBs->end(); iBl++)
                        {
                            if (iBl->find("_symmetry_space_group_") !=std::string::npos)
                            {
                                getCifCrystInfo(iBs);
                                break;
                            }
                            else if (iBl->find("_symmetry_equiv_pos_as_xyz") !=std::string::npos)
                            {
                                getCifSymOps(iBs);
                                break;
                            }
                            else if (iBl->find("_atom_site_label") !=std::string::npos
                                     && iBl->size()==16)
                            {
                                // std::cout << "_atom_site_label" << std::endl;
                                getCifAtomInfo(iBs);
                                break;
                            }
                        }
               
                    }
                
                }
            }
            if (itsCurCryst !=NullPoint)
            {
                allCryst.push_back(*itsCurCryst);
            }
                 
            // check
           
            std::cout << "There are " << (int)allAtoms.size() 
                    << " atoms in the system. They are " << std::endl;
            
            for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                   iA !=allAtoms.end(); iA++)
            {
                std::cout << "Atom " << iA->id << std::endl
                          << "Its type symbol " << iA->chemType << std::endl
                          << "Its coordinates (fractional) " << std::endl
                          << "x:   " << iA->fracCoords[0] << std::endl
                          << "y:   " << iA->fracCoords[1] << std::endl
                          << "z:   " << iA->fracCoords[2] << std::endl
                          << "Its coordinates  " << std::endl
                          << "x:   " << iA->coords[0] << std::endl
                          << "y:   " << iA->coords[1] << std::endl
                          << "z:   " << iA->coords[2] << std::endl;           
            }
        }
        
    } 
    
    void GenCifFile::setupSystem2()
    {
       if (inFile.is_open() )
        { 
            std::vector<std::vector<std::string> > tBlocs;
            // make sure
            itsCurBlock = "";
            std::string tRecord="";
            
            std::vector<std::string>     tAllLines;
            std::vector<std::string>     tRAllLines;
            std::vector<std::string>     tBlocLines;
            
            
            while(!inFile.eof())
            {   
                std::getline(inFile, tRecord);
                tRecord = TrimSpaces(tRecord);
                tAllLines.push_back(tRecord);
                // std::cout <<  tRecord << std::endl;
                // only use a few blocks in the cif file
                if (checkR)
                {
                    if (tRecord.find("_refine_ls_R_factor_all") !=std::string::npos)
                    {
                        std::cout << "R factor line: " << tRecord << std::endl;
                        tRAllLines.push_back(tRecord);
                    }
                
                    if (tRecord.find("_refine_ls_R_factor_gt") !=std::string::npos)
                    {
                        std::cout << "R factor line: " << tRecord << std::endl;
                        tRAllLines.push_back(tRecord);
                    }
                
                    if (tRecord.find("_refine_ls_R_factor_obs") !=std::string::npos)
                    {
                        std::cout << "R factor line: " << tRecord << std::endl;
                        tRAllLines.push_back(tRecord);
                    }
                }
                
                if (tRecord.find("loop_") !=std::string::npos 
                    || (tRecord.find("data_") !=std::string::npos
                        && tRecord.find("_data_") ==std::string::npos))
                {
                    if (!tBlocLines.empty())
                    {
                        tBlocs.push_back(tBlocLines);
                        tBlocLines.clear();
                    }
                    
                    itsCurBlock = "loop";
                }
                else if (itsCurBlock=="loop" && tRecord.size() >0 )
                {
                    if (tRecord[0] !='#')
                    {
                        tBlocLines.push_back(tRecord);
                    }
                }
            }
            
            if (!tBlocLines.empty())
            {
                tBlocs.push_back(tBlocLines);
                tBlocLines.clear();
            }
            
            inFile.close();
            if (checkR)
            {
                std::cout << "tRAllLine.size " <<  tRAllLines.size() << std::endl;
                checkRFact(tRAllLines);
            }
            else
            {
                RFactorOK = true;
            }
            
            
            
            checkPowder(tAllLines);
            
            checkNeutronD(tAllLines);
            
            if (!notPowder)
            {
                return ;
            }
           
            
            // std::cout << "Number of data blocks in the cif file is " 
            // << tBlocs.size() << std::endl;
            
            /*
            for (std::vector<std::vector<std::string> >::iterator 
                 iBloc=tBlocs.begin();
                 iBloc !=tBlocs.end(); iBloc++)
            {
                std::cout << "====One block lines " << std::endl;
                for (std::vector<std::string>::iterator iLine=iBloc->begin();
                        iLine != iBloc->end(); iLine++)
                {
                    std::cout << *iLine << std::endl;
                }
            }
            */
            
            // 
            //if ((int)tBlocs.size()!=0 && RFactorOK )
            if ((int)tBlocs.size()!=0)
            {
                
                std::map<std::string, std::string>   rowProps;
                std::map<int, std::map<ID, std::vector<std::string> > > colProps;
           
                int idxB=0;
                for (std::vector<std::vector<std::string> >::iterator iBs=
                        tBlocs.begin(); iBs !=tBlocs.end(); iBs++)
                {
                    
                    if((int)iBs->size() !=0)
                    {
                        getPropsToMaps(iBs, rowProps, colProps, idxB);
                    }
                }
                 
                
                // Check 
                
                /*
                std::cout << "There are " << rowProps.size() 
                          << " properties. They are: " << std::endl;
                for (std::map<std::string, std::string>::iterator iRPs=rowProps.begin();
                        iRPs !=rowProps.end(); iRPs++)
                {
                    std::cout << iRPs->first << " -> " << iRPs->second << std::endl;
                }
                
                
                std::cout << std::endl << "There are " << colProps.size() 
                          << " blocks of properties. They are: " << std::endl;
                
                for (std::map<int, std::map<ID, std::vector<std::string> > >::iterator 
                        iCPs=colProps.begin(); iCPs != colProps.end(); iCPs++)
                {
                    std::cout << "For block " << iCPs->first << " : " << std::endl;
                    std::cout << "Property labels: " << std::endl;
                    for (std::vector<std::string>::iterator 
                          iOneB=iCPs->second["lab"].begin(); 
                          iOneB != iCPs->second["lab"].end(); iOneB++)
                    {
                        std::cout << (*iOneB) << std::endl;
                    }
                    std::cout << "The lines associated : " << std::endl;
                    for (std::vector<std::string>::iterator 
                          iOneB=iCPs->second["cont"].begin(); 
                          iOneB != iCPs->second["cont"].end(); iOneB++)
                    {
                        std::cout << (*iOneB) << std::endl;
                    }
                }
                */
               
                
                // Now select what we need from the input cif file.
                selectPropsToMaps(rowProps, colProps);
                
                if (!lErr)
                {
                    checkAtomElementID();
                
                    getHAtomIdxs();
                    
                    std::cout << "there are " << allHAtomIdx.size() 
                              << "H atom in the initial atom set " << std::endl;
                    for (std::vector<int>::iterator iH=allHAtomIdx.begin();
                            iH != allHAtomIdx.end(); iH++)
                    {
                        std::cout << "Atom " << allAtoms[*iH].id 
                                  << " is a H atom" << std::endl;
                    }
                    
                   
                    
                    // check if the assigned atom element types are right ones
                    PeriodicTable aPTab;
                    for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                        iA !=allAtoms.end(); iA++)
                    {
                        if (!assignElementType(aPTab, iA->chemType, iA))
                        {
                            std::cout << "Warn: atom " << iA->id 
                                      << " has element type "
                                      << iA->chemType 
                                      << " which is not in the Periodic Table"
                                      << std::endl;
                        }
                    }
                
                    
                    checkNonHAtomOccp();
                   
                    setAtomsMetalType();
                    setAtomOxiType();
                    
                    checkCalcAtoms();
                
                     // Check
                
                    std::cout << "There are " << (int)allAtoms.size() 
                              << " atoms in the system. They are " << std::endl;
            
                    for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                          iA !=allAtoms.end(); iA++)
                    {
                        std::cout << "Atom " << iA->id << std::endl
                                  << " serial Number " << iA->seriNum << std::endl
                                  << "Its type symbol " << iA->chemType << std::endl
                                  << "Is it metal ";
                   
                        if (iA->isMetal)
                        {
                            std::cout << " Yes " << std::endl;
                        }
                        else
                        {
                            std::cout << " No " << std::endl;
                        }
                        
                        std::cout << "Is that atom from calculations ? ";
                        if (iA->fromCalc)
                        {
                            std::cout <<  "Yes" << std::endl;
                        }
                        else
                        {
                            std::cout << "No, it is determined from experimental data "
                                      << std::endl;
                        }
                        
                        std::cout << "Its form charge " << iA->formalCharge << std::endl;
                        std::cout << "Its occupancy " << iA->ocp << std::endl;
                        std::cout << "Its coordinates (fractional) " << std::endl
                                  << "x:   " << iA->fracCoords[0] << std::endl
                                  << "y:   " << iA->fracCoords[1] << std::endl
                                  << "z:   " << iA->fracCoords[2] << std::endl;
                                  //<< "Its coordinates  " << std::endl
                                  //<< "x:   " << iA->coords[0] << std::endl
                                  //<< "y:   " << iA->coords[1] << std::endl
                                  //<< "z:   " << iA->coords[2] << std::endl;  
                   
                        //std::cout << "Is that atom in read-in unit ? " 
                        //          << iA->isInPreCell << std::endl;
                   
                    }
                }
            }
       }
      
    }
    
    
    void GenCifFile::setupSystemCSD()
    {
       if (inFile.is_open() )
        { 
            std::vector<std::vector<std::string> > tBlocs;
            std::map<std::string, std::string>     all2Cols;
            
            // make sure
            itsCurBlock = "";
            std::string tRecord="";
            
            std::vector<std::string>     tAllLines;
            std::vector<std::string>     tBlocLines;
            REAL tROK1 = -1.0;
            REAL tROK2 = -1.0;
            
            while(!inFile.eof())
            {   
                std::getline(inFile, tRecord);
                tRecord = TrimSpaces(tRecord);
                std::vector<std::string> aPair;
                StrTokenizeGen(tRecord, aPair);
                if (aPair.size()==2)
                {
                    std::string tP0 = TrimSpaces(aPair[0]);
                    std::string tP1 = TrimSpaces(aPair[1]);
                    if (tP0.size() > 0)
                    {
                        if (tP0[0] =='_')
                        {
                            all2Cols[tP0] = tP1;
                        }
                        else if (tP0[0] !='#')
                        {
                            tAllLines.push_back(tRecord);
                        }
                    }
                }
                else if (aPair.size() > 0)
                {
                    tAllLines.push_back(tRecord);
                }
            }
            
            for (unsigned i=0; i < tAllLines.size(); i++)
            {   
                tRecord = tAllLines[i];
                tRecord = TrimSpaces(tRecord);
                // std::cout <<  tRecord << std::endl;
                std::vector<std::string> tBuf;
                std::vector<std::string> tBuf_t;
                
                // only use a few blocks in the cif file
                
                if (tRecord.find("loop_") !=std::string::npos 
                    || (tRecord.find("data_") !=std::string::npos
                        && tRecord.find("_data_") ==std::string::npos))
                {
                    if (!tBlocLines.empty())
                    {
                        tBlocs.push_back(tBlocLines);
                        tBlocLines.clear();
                    }
                    
                    itsCurBlock = "loop";
                }
                else if (itsCurBlock=="loop" && tRecord.size() >0 )
                {
                    if (tRecord[0] !='#')
                    {
                        tBlocLines.push_back(tRecord);
                    }
                }
            }
            
            if (!tBlocLines.empty())
            {
                tBlocs.push_back(tBlocLines);
                tBlocLines.clear();
            }
            
            if (tROK1 > 0.0 && tROK2 > 0.0)
            {
                if (tROK1 <RTHRESHOLD_U && tROK2 < RTHRESHOLD_U)
                {
                    RFactorOK = true;
                }
            }
            else if (tROK1 > 0.0 && tROK2 < 0.0)
            {
                if (tROK1 <RTHRESHOLD_U)
                {
                    RFactorOK = true;
                }
            }
            else if (tROK1 < 0.0 && tROK2 > 0.0)
            {
                if (tROK2 <RTHRESHOLD_U)
                {
                    RFactorOK = true;
                }
            }
            
            if (!RFactorOK)
            {
                errMsg.push_back("R FACTOR IS TOO HIGH\n");
            }
            
            inFile.close();
           
            checkPowder(tAllLines);
            checkNeutronD(tAllLines);
            
            if (!notPowder)
            {
                return ;
            }
            // std::cout << "Number of data blocks in the cif file is "
            // << tBlocs.size() << std::endl;
            
            /*
            for (std::vector<std::vector<std::string> >::iterator 
                    iBloc=tBlocs.begin();
                    iBloc !=tBlocs.end(); iBloc++)
            {
                std::cout << "====One block lines " << std::endl;
                for (std::vector<std::string>::iterator iLine=iBloc->begin();
                        iLine != iBloc->end(); iLine++)
                {
                    std::cout << *iLine << std::endl;
                }
            }
            */
            
            // 
            //if ((int)tBlocs.size()!=0 && RFactorOK )
            if ((int)tBlocs.size()!=0)
            {
                
                std::map<std::string, std::string>   rowProps;
                std::map<int, std::map<ID, std::vector<std::string> > > colProps;
           
                int idxB=0;
                for (std::vector<std::vector<std::string> >::iterator iBs=
                        tBlocs.begin(); iBs !=tBlocs.end(); iBs++)
                {
                    
                    if((int)iBs->size() !=0)
                    {
                        getPropsToMaps(iBs, rowProps, colProps, idxB);
                    }
                }
                 
                
                // Check 
                
                /*
                std::cout << "There are " << rowProps.size() 
                          << " properties. They are: " << std::endl;
                for (std::map<std::string, std::string>::iterator iRPs=rowProps.begin();
                        iRPs !=rowProps.end(); iRPs++)
                {
                    std::cout << iRPs->first << " -> " << iRPs->second << std::endl;
                }
                
                
                std::cout << std::endl << "There are " << colProps.size() 
                          << " blocks of properties. They are: " << std::endl;
                
                for (std::map<int, std::map<ID, std::vector<std::string> > >::iterator 
                        iCPs=colProps.begin(); iCPs != colProps.end(); iCPs++)
                {
                    std::cout << "For block " << iCPs->first << " : " << std::endl;
                    std::cout << "Property labels: " << std::endl;
                    for (std::vector<std::string>::iterator 
                          iOneB=iCPs->second["lab"].begin(); 
                          iOneB != iCPs->second["lab"].end(); iOneB++)
                    {
                        std::cout << (*iOneB) << std::endl;
                    }
                    std::cout << "The lines associated : " << std::endl;
                    for (std::vector<std::string>::iterator 
                          iOneB=iCPs->second["cont"].begin(); 
                          iOneB != iCPs->second["cont"].end(); iOneB++)
                    {
                        std::cout << (*iOneB) << std::endl;
                    }
                }
                */
               
               
                // Now select what we need from the input cif file.
                selectPropsToMaps(rowProps, colProps);
                
                if (!lErr)
                {
                    checkAtomElementID();
                
                
                    // check if the assigned atom element types are right ones
                    PeriodicTable aPTab;
                    for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                        iA !=allAtoms.end(); iA++)
                    {
                        if (!assignElementType(aPTab, iA->chemType, iA))
                        {
                            std::cout << "Warn: atom " << iA->id << " has element type "
                                      << iA->chemType << " which is not in the Periodic Table"
                                      << std::endl;
                        }
                    }
                
                    checkNonHAtomOccp();
                    setAtomsMetalType();
                    setAtomOxiType();
                    
                    checkCalcAtoms();
                
                     // Check
                
                    std::cout << "There are " << (int)allAtoms.size() 
                              << " atoms in the assym unit. They are " << std::endl;
            
                    for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                          iA !=allAtoms.end(); iA++)
                    {
                        std::cout << "Atom " << iA->id << std::endl
                                  << " serial Number " << iA->seriNum << std::endl
                                  << "Its type symbol " << iA->chemType << std::endl
                                  << "Is it metal ";
                   
                        if (iA->isMetal)
                        {
                            std::cout << " Yes " << std::endl;
                        }
                        else
                        {
                            std::cout << " No " << std::endl;
                        }
                        
                        /*
                        std::cout << "Is that atom from calculations ? ";
                        if (iA->fromCalc)
                        {
                            std::cout <<  "Yes" << std::endl;
                        }
                        else
                        {
                            std::cout << "No, it is determined from experimental data "
                                      << std::endl;
                        }
                        */
                        
                        std::cout << "Its form charge " << iA->formalCharge << std::endl;
                        std::cout << "Its occupancy " << iA->ocp << std::endl;
                        std::cout << "Its coordinates (fractional) " << std::endl
                                  << "x:   " << iA->fracCoords[0] << std::endl
                                  << "y:   " << iA->fracCoords[1] << std::endl
                                  << "z:   " << iA->fracCoords[2] << std::endl;
                                  //<< "Its coordinates  " << std::endl
                                  //<< "x:   " << iA->coords[0] << std::endl
                                  //<< "y:   " << iA->coords[1] << std::endl
                                  //<< "z:   " << iA->coords[2] << std::endl;  
                   
                        //std::cout << "Is that atom in read-in unit ? " 
                        //          << iA->isInPreCell << std::endl;
                   
                    }
                }
            }
       }
       
    }
    
    
    void GenCifFile::checkRFact(std::vector<std::string>& tLines)
    {
        
        REAL tROK1 = -1.0;
        REAL tROK2 = -1.0;
        REAL tROK3 = -1.0;
        
        for (unsigned i=0; i < tLines.size(); i++)
        {
            std::vector<std::string> tBuf; 
            if (tLines[i].find("_refine_ls_R_factor_all") !=std::string::npos)
            {    
                std::cout << "R factor line: " << tLines[i] << std::endl;
                            
                StrTokenize(tLines[i],  tBuf);
                    
                if (tBuf.size()==2)
                {
                    tROK1=StrToReal(tBuf[1]);
                    if (tROK1 >=1.0)
                    {
                        tROK1 = tROK1/100.0;
                    }
                    std::cout << "_refine_ls_R_factor_all="
                              << tROK1 << std::endl;
                }
            }
            else if (tLines[i].find("_refine_ls_R_factor_gt") 
                     !=std::string::npos)
            {
                std::cout << "R factor line: " << tLines[i] << std::endl;
                StrTokenize(tLines[i],  tBuf);
                if (tBuf.size()==2)
                {
                    tROK2=StrToReal(tBuf[1]);
                    if (tROK2 >=1.0)
                    {
                        tROK2 = tROK2/100.0;
                    }
                    std::cout << "_refine_ls_R_factor_gt="
                              << tROK2 << std::endl;
                }
            }
            else if (tLines[i].find("_refine_ls_R_factor_obs") 
                     !=std::string::npos)
            {
                std::cout << "R factor line: " << tLines[i] << std::endl;
                    
                StrTokenize(tLines[i],  tBuf);
                if (tBuf.size()==2)
                {
                    tROK3=StrToReal(tBuf[1]);
                    if (tROK3 >=1.0)
                    {
                        tROK3 = tROK3/100.0;
                    }
                    std::cout << "_refine_ls_R_factor_all="
                              << tROK3 << std::endl;
          
                }
            }
        }

        if ((tROK1 > 0.0 && tROK2 > 0.0) || tROK3 >0.0)
        {
                
            if (tROK1 <RTHRESHOLD_U && 
                (tROK2 < RTHRESHOLD_U || tROK3 < RTHRESHOLD_U))
            {
                    RFactorOK = true;                
            }
        }
        else if (tROK1 > 0.0 && tROK2 < 0.0)
        {
            if (tROK1 <RTHRESHOLD_U)
            {
                    RFactorOK = true;
            }
        }
        else if (tROK1 < 0.0 && tROK2 > 0.0)
        {
            if (tROK2 <RTHRESHOLD_U)
            {
                RFactorOK = true;
            }
        }
            
        if (!RFactorOK)
        {
            lErr = true;
            errMsg.push_back("R FACTOR IS TWO HIGH\n");
        }
        
    }
    
    void GenCifFile::checkPowder(std::vector<std::string>& tLines)
    {
        if(tLines.size() !=0)
        {
            for (std::vector<std::string>::iterator iL=tLines.begin();
                    iL != tLines.end(); iL++)
            {
                StrLower(*iL);
                // std::cout << *iL << std::endl;
                if (iL->find("#")==iL->npos && (iL->find("powder") !=iL->npos
                       || iL->find("_ph_") !=iL->npos))
                {
                    notPowder = false;
                    errMsg.push_back("POWDER DIFFRACTIONS DETECTED\n");
                    break;
                }
            }
        }  
    }
    
    void GenCifFile::checkNeutronD(std::vector<std::string>& tLines)
    {
        if(tLines.size() !=0)
        {
            for (std::vector<std::string>::iterator iL=tLines.begin();
                    iL != tLines.end(); iL++)
            {
                StrLower(*iL);
                // std::cout << *iL << std::endl;
                if (iL->find("#")==iL->npos && (iL->find("neutron") !=iL->npos))
                {
                    notNeuD = false;
                    errMsg.push_back("NEUTRON DIFFRACTIONS DETECTED\n");
                    break;
                }
            }
        }  
    }
    
    void GenCifFile::checkCalcAtoms()
    {
        if (allAtoms.size() >0)
        {
            std::vector<int> idxsCAt;
            for (std::vector<AtomDict>::iterator iAt=allAtoms.begin();
                iAt !=allAtoms.end(); iAt++)
            {
                if (iAt->fromCalc && iAt->chemType !="H" && iAt->chemType !="D")
                {
                    hasHeavyCalcAtoms = true;
                    //std::cout << iAt->id << " iAt->fromCalc " 
                    //          << iAt->fromCalc << std::endl;
                    idxsCAt.push_back(iAt->seriNum);
                }
            }
            
            //std::cout << "The system hasHeavyCalcAtoms? "
            //          << eavyCalcAtoms << std::endl;
            if (idxsCAt.size() !=0)
            {
                for (unsigned i=0; i < idxsCAt.size(); i++)
                {
                    errMsg.push_back("Reject structures : \n");
                    errMsg.push_back("Atom " + allAtoms[idxsCAt[i]].id 
                                     + " is from calculations\n" );
                }
            }
        } 
    }
    
    
    
    bool GenCifFile::checkOverAll(int tWorkMode)
    {
        bool aRet= true;
        std::cout << "not neutron diffr " << notNeuD << std::endl;
        std::cout << "not powder diffr " << notPowder << std::endl;
        std::cout << "colidOK " << colidOK << std::endl;
        std::cout << "hasHeavyCalcAtoms " << hasHeavyCalcAtoms << std::endl;
        std::cout << "RFactorOK " << RFactorOK << std::endl;
        
        if (tWorkMode==3111)
        {
            
            if (!notPowder || notNeuD
                || !colidOK || hasHeavyCalcAtoms ||!symmOprOK|| lErr)
            {
                aRet=false;
            }
        }
        else if (tWorkMode==31 
             || tWorkMode==311
             || tWorkMode ==312
             || tWorkMode ==313 
             || tWorkMode ==314
             || tWorkMode==32 
             || tWorkMode==33)
        {
            if (!notPowder || !resolOK
                || !RFactorOK || !colidOK 
                || hasHeavyCalcAtoms ||!symmOprOK|| lErr)
            {
                aRet=false;
            }
        }
        //std::cout << aRet << std::endl;
        return aRet;
    }
    
    
    void GenCifFile::setAllCrits(std::vector<std::string>& tAllLines)
    {
        for (std::vector<std::string>::iterator iL=tAllLines.begin();
               iL !=tAllLines.end(); iL++)
        {
            std::cout << *iL << std::endl;
            if (iL->size() > 0)
            {
                std::vector<std::string> tBuf;
                StrTokenize(*iL,  tBuf, ':');
                if(tBuf.size()==2)
                {
                    StrUpper(tBuf[0]);
                    if (tBuf[0].find("NOCHECK") !=tBuf[0].npos)
                    {
                        StrUpper(tBuf[1]);
                        if (tBuf[1].find("YES") !=tBuf[1].npos)
                        {
                          nonCheck = true;
                        }
                    }
                    else if(tBuf[0].find("CHECKRFACT") !=tBuf[0].npos)
                    {
                        RTHRESHOLD_U=StrToReal(tBuf[1]);
                        std::cout << RTHRESHOLD_U << std::endl;
                        if (RTHRESHOLD_U < 0.0)
                        {
                            checkR = false;
                        }
                    }
                    else if (tBuf[0].find("CHECKRESOL") !=tBuf[0].npos)
                    {
                        RESOLTHRESHOLD_U=StrToReal(tBuf[1]);
                        if (RESOLTHRESHOLD_U < 0.0)
                        {
                            checkResol = false;
                        }
                        std::cout << RESOLTHRESHOLD_U << std::endl;
                    }
                }
            }
        }
        
        if (nonCheck)
        {
            checkR     = false;
            checkResol = false;
        }
        std::cout << "check factor ? " << checkR << std::endl;
        std::cout << "check resol ? " << checkResol << std::endl;
        
        
    }
    
    void GenCifFile::checkAtomElementID()
    {
        PeriodicTable aPTab;
       
        for (std::vector<AtomDict>::iterator iAt=allAtoms.begin();
                iAt !=allAtoms.end(); iAt++)
        {
            
            //std::cout << "Atom id " << iAt->id << std::endl;
            //std::cout << "Atom element symbol " << iAt->chemType << std::endl;
            /*
            if (iAt->chemType.empty())
            {
                std::cout << "empty " << iAt->chemType.empty() << std::endl;
            }
            
            if (aPTab.elements.find(iAt->chemType) == aPTab.elements.end())
            {
                std::cout << "Not find " << std::endl;
            }
            */
            
            if (iAt->chemType.empty() || ((iAt->chemType.size() < 4) &&
                aPTab.elements.find(iAt->chemType) == aPTab.elements.end()))
                     
            {
                bool tFind = false;
                
                if (!iAt->id.empty())
                {
                    if (iAt->id.size() >=2)
                    {
                        if(aPTab.elements.find(iAt->id.substr(0,2)) 
                                !=aPTab.elements.end())
                        {
                            iAt->chemType=iAt->id.substr(0,2);
                            tFind = true;
                        }
                        else if (aPTab.elements.find(iAt->id.substr(0,1))
                                  !=aPTab.elements.end())
                        {
                            iAt->chemType=iAt->id.substr(0,1);
                            tFind = true;
                        }
                    }
                    else if(iAt->id.size()==1)
                    {
                        if (aPTab.elements.find(iAt->id)
                              !=aPTab.elements.end())
                        {
                            iAt->chemType=iAt->id;
                            tFind = true;
                        }
                        else if(iAt->id.compare(".")==0)
                        {
                            iAt->chemType = "Po";
                            tFind  = true;
                        }
                    }
                    
                }
                
                if (!tFind)
                {
                    std::cout << "Atom " << iAt->id << "'s element type could not be handled with in this program "
                              << std::endl;
                    std::string aSL = "Atom " + iAt->id + "'s element type could not be handled with in this program \n";
                    errMsg.push_back(aSL);
                }
            }
            else if (iAt->chemType.size() >2)
            {        
                if (iAt->chemType.find("+") !=std::string::npos
                    || iAt->chemType.find("-") !=std::string::npos)
                {
                    std::string aSC, aSN;
                    for (unsigned i =0; i < iAt->chemType.size(); i++)
                    {
                        if(std::isdigit(iAt->chemType[i]))
                        {
                           aSN+=(iAt->chemType[i]); 
                        }
                        else if (iAt->chemType[i]=='-')
                        {
                            aSN = iAt->chemType[i]+ aSN;
                        }
                        else if (iAt->chemType[i] !='+')
                        {
                            aSC+=(iAt->chemType[i]);
                        }
                    }
                    
                    
                    
                    if (aSC.size() < 3)
                    {
                        iAt->chemType=aSC;
                    }
                    else
                    {
                        std::cout << "Bug: atom " << iAt->id << " is a element "
                                  << iAt->chemType << std::endl;
                    }
                    
                    if(aSN.size() !=0)
                    {
                        iAt->formalCharge=StrToReal(aSN);
                    }
                    
                }
            }

            // std::cout << "Its atom element type " << iAt->chemType << std::endl;
            
        }
    }
    
    
    void GenCifFile::checkNonHAtomOccp()
    {
        float nAll = 0.0, nColid=0.0;
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                     iA !=allAtoms.end(); iA++)
        {
            if (iA->chemType.find("H")==std::string::npos)
            {
                nAll +=1.0;
                
                if (iA->ocp < 0.99)
                {
                    nColid+=1.0;
                }
            }
        }
        
        //if (nColid/nAll >0.2)
        if (nColid > 0)
        {
            //std::cout << "nAll " << nAll << std::endl;
            //std::cout << "nColid " << nColid << std::endl;
            //std::cout << "ratio " << nColid/nAll << std::endl;
            //errMsg.push_back("HIGHLY COLLIDED STRUCUTRE \n");
            
            
            errMsg.push_back("Some non-H atoms have less than 1.0 occupancy \n");
            colidOK = false;
        }
        
    }
    
    void GenCifFile::setAtomsMetalType()
    {
        std::vector<ID>  allMetals;
        std::vector<ID>  pureMetals;
        std::vector<ID>  allMetalloids;
        
        initMetalTab(pureMetals);
        initMetalloidTab(allMetalloids);
        for (std::vector<ID>::iterator iM=pureMetals.begin();
                iM !=pureMetals.end(); iM++)
        {
            allMetals.push_back(*iM);
        }
        for (std::vector<ID>::iterator iM=allMetalloids.begin();
                iM !=allMetalloids.end(); iM++)
        {
            allMetals.push_back(*iM);
        }
        
        //std::cout << "Number of metal in the table " << allMetals.size() 
        //          << std::endl;
        //std::cout << "Number of metalloids in the table " 
        //          << allMetalloids.size() << std::endl;
        
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                     iA !=allAtoms.end(); iA++)
        {
            if (isMetal(pureMetals, iA->chemType))
            {
                iA->matType = 1;
            }
            else if (isMetalloid(allMetalloids, iA->chemType))
            {
                iA->matType = 2;
            }
            iA->isMetal = isMetal(allMetals, iA->chemType);
            if (iA->isMetal)
            {
                if (!hasMetal)
                {
                    hasMetal = true;
                }
            }
            /*
            std::string TP="organic ";
            
            if (iA->matType ==1)
            {
                TP = "metal";
            }
            else if (iA->matType == 2)
            {
                TP = "metalloid";
            }
            
            std::cout << "Atom " << iA->id << " of " << iA->chemType
                      << " is a " << TP << " atom " << std::endl;
             */
        }
        
        
    }
    
    void GenCifFile::getPropsToMaps(std::vector<std::vector<std::string> >::iterator tOneBlockLines, 
                                    std::map<std::string,std::string>       & tRowProps, 
                                    std::map<int,std::map<ID,std::vector<std::string> > > & tColProps,
                                    int & tIdxB)
    {
        bool lLab=false;
        // std::cout << "A block has been read " << std::endl;
        for (std::vector<std::string>::iterator iBl=tOneBlockLines->begin();
                               iBl != tOneBlockLines->end(); iBl++)
        {
            (*iBl) = TrimSpaces((*iBl));
            //std::cout << "The lines is : " 
            //          << (*iBl) << std::endl;
            std::vector<std::string> tBuf0, tBuf;
            
            if (iBl->find('\'') != std::string::npos && iBl->substr(0,1).compare("_")==0)
            {
                StrTokenize(*iBl, tBuf0, '\'');
                // std::cout << "line " << (*iBl) << std::endl;
                // std::cout << "strgrp size " << tBuf0.size() << std::endl;
                for (unsigned i0=0; i0 < tBuf0.size(); i0++)
                {
                    // std::cout << "string " << i0+1 << tBuf0[i0] << std::endl;
                    if (!tBuf0[i0].empty())
                    {
                        tBuf.push_back(tBuf0[i0]);
                    }
                }
                /*
                for (unsigned i0=0; i0 < tBuf.size(); i0++)
                {
                    std::cout << "after string  " << i0+1 << tBuf[i0] << std::endl;
                }
                 */
            }
            else if (iBl->find('\"') != std::string::npos)
            {
                
                std::vector<std::string> tBuf0, tBuf00;
                StrTokenize(*iBl, tBuf00);
                if (tBuf00[0].find('\"') != std::string::npos)
                {
                    for (unsigned i0=0; i0 < tBuf00.size(); i0++)
                    {
                        if (tBuf00[i0].size() >0)
                        {
                            tBuf.push_back(tBuf00[i0]);
                        }
                    }
                }
                else
                {
                    StrTokenize(*iBl, tBuf0, '\"');
                    for (unsigned i0=0; i0 < tBuf0.size(); i0++)
                    {
                        if (tBuf0[i0].size() >0)
                        {
                            tBuf.push_back(tBuf0[i0]);
                        }
                    }
                }
            }
            else if (iBl->find('\'') != std::string::npos)
            {
                std::vector<std::string> tBuf0, tBuf00;
                StrTokenize(*iBl, tBuf00);
                if (tBuf00[0].find('\'') != std::string::npos)
                {
                    for (unsigned i0=0; i0 < tBuf00.size(); i0++)
                    {
                        if (tBuf00[i0].size() >0)
                        {
                            tBuf.push_back(tBuf00[i0]);
                        }
                    }
                }
                else
                {
                    StrTokenize(*iBl, tBuf0, '\'');
                    for (unsigned i0=0; i0 < tBuf0.size(); i0++)
                    {
                        if (tBuf0[i0].size() >0)
                        {
                            tBuf.push_back(tBuf0[i0]);
                        }
                    }
                }
                
                // tBuf.push_back((*iBl));
            }
            else
            {
                
                StrTokenize(*iBl, tBuf);
            }
            
            
            if (tBuf.size()==1)
            {
                if (tBuf[0].find("_") !=std::string::npos)
                { 
                    
                    if (lLab) 
                    {
                        tColProps[tIdxB]["lab"].push_back(TrimSpaces(tBuf[0]));
                        
                    }
                    else
                    {
                        
                        lLab=true;
                        tIdxB++;
                        tColProps[tIdxB]["lab"].push_back(TrimSpaces(tBuf[0]));
                        
                    }
                }
                else
                {
                    if (lLab)
                    {
                        lLab  = false;
                    }
                    tColProps[tIdxB]["cont"].push_back((*iBl));   
                    
                }
                
            }
            else if (tBuf.size() ==2 && tBuf[0].find("_") !=std::string::npos)
            {
                tRowProps[TrimSpaces(tBuf[0])] = TrimSpaces(tBuf[1]);
                
            }
            else if (tBuf.size() > 1 && tBuf[0][0] !='_')
            {
                // if(tBuf.size()==tColProps[tIdxB]["lab"].size())
                //{
                    if (lLab)
                    {
                        lLab  = false;
                    }
                    //std::cout << "the firs label is " << tColProps[tIdxB]["lab"][0]
                    //           << std::endl;
                    tColProps[tIdxB]["cont"].push_back((*iBl));
                    //std::cout << "cont line: " << *iBl << std::endl;
                // }
            }
        } 
    }
    
    void GenCifFile::selectPropsToMaps(std::map<std::string,std::string> 
                                       & tRowProps, 
                                       std::map<int,std::map<ID,
                                       std::vector<std::string> > > & tColProps)
    {   
        
        getCifCrystInfo(tRowProps, tColProps);
      
        bool lSymOps=false, lAtms=false, lTypes=false;
        for (std::map<int,std::map<ID,std::vector<std::string> > >::iterator 
              iBl=tColProps.begin(); iBl != tColProps.end(); iBl++)
        {
            if (!lSymOps)
            {
                /*
                for (std::vector<std::string>::iterator 
                     iM=iBl->second["lab"].begin(); 
                     iM != iBl->second["lab"].end(); iM++)
                {
                    std::cout << "print lab " << *iM << std::endl;
                }
                */
                int aPos=getKeyWordPos("_symmetry_equiv_pos_as_xyz", 
                                            iBl->second["lab"]);
                if (aPos ==-1)
                {
                                      
                    aPos=getKeyWordPos("_space_group_symop_operation_xyz",
                                            iBl->second["lab"]);
                }
                if (aPos !=-1)
                {
                    lSymOps=true;
                    getCifSymOps(iBl->second);
                }
               
            }
            
            
            if(!lAtms)
            {
                
                if (getKeyWordPos("_atom_site_type_symbol", iBl->second["lab"]) !=-1
                    || getKeyWordPos("_atom_site_label", iBl->second["lab"]) !=-1)
                {
                    
                    lAtms = true;
                    getCifAtomInfo(iBl->second);
                }
            }
          
            
            if(!lTypes)
            {
                if (getKeyWordPos("_atom_type_symbol", iBl->second["lab"]) !=-1)
                {
                    lTypes = true;
                    getCifAtomOxiInfo(iBl->second);
                }
            }
        }
        
        if (itsCurCryst !=NullPoint)
        {
            
            allCryst.push_back(*itsCurCryst);
            // std::cout << "number of crystal : " << allCryst.size() << std::endl;
        }
        else
        {
            std::cout << "No crystal " << std::endl;
            errMsg.push_back("No crystal ");
            lErr = true;   
        }
        
        if (!lSymOps)
        {
            symmOprOK = false;
            std::cout << "No symmetry operators, can not generator molecules "
                      << std::endl;
            errMsg.push_back("No symmetry operators, can not generator molecules \n");
            lErr = true;
        }
    }
    
    void GenCifFile::initAllCifKeys()
    {
        //std::string clibMonDir(std::getenv("CLIBD_MON"));
        //std::string fName(clibMonDir);
        //fName.append("/list/cif_tag.list");
        std::string clibMonDir(std::getenv(""));
        std::string fName(clibMonDir);
        fName.append("/tables/cif_tag.list");
        
        std::ifstream fCifKeys(fName.c_str());
        if (fCifKeys.is_open())
        {
            
            while (!fCifKeys.eof())
            {
                std::string tRecord="";
                std::getline(fCifKeys, tRecord);
                TrimSpaces(tRecord);
                existCifKeys[tRecord]=false;
            }
            fCifKeys.close();
        }
    }
    
    void GenCifFile::getCifCrystInfo(std::vector<std::vector<std::string> >::iterator tBs)
    {
        // std::cout << "Get cell parameters and space group " << std::endl;
        itsCurCryst = new CrystInfo();
        itsCurCryst->itsSpaceGroup = new SpaceGroupMember();
        itsCurCryst->itsCell       = new Cell();
        
        for (std::vector<std::string>::iterator iBl=tBs->begin();
                iBl!=tBs->end(); iBl++)
        {   
            if (iBl->find("_space_group_IT_number") !=std::string::npos)
            {
                std::vector<std::string> tBuf;
                StrTokenize(*iBl, tBuf);
                if ((int)tBuf.size() == 2)
                {
                    itsCurCryst->itsSpaceGroup->sgNum = StrToInt(tBuf[1]); 
                }
                else
                {
                    std::cout << "No space group number in line " << *iBl << std::endl;
                    lErr = true;
                    errMsg.push_back("No space group number in line");
                    exit(1);
                }
            }
            else if (iBl->find("_symmetry_cell_setting") !=std::string::npos)
            {
                std::vector<std::string> tBuf;
                StrTokenize(*iBl, tBuf);
                if ((int)tBuf.size() == 2)
                {
                    itsCurCryst->itsCell->lattice = tBuf[1]; 
                }
                else
                {
                    std::cout << "Wrong Line format for the lattice in line" << *iBl << std::endl;
                    errMsg.push_back("Wrong Line format for the lattice in line");
                    lErr=true;
                }
            }
            else if (iBl->find("_symmetry_space_group_name_Hall") !=std::string::npos)
            {
                std::vector<std::string> tBuf;
                if (iBl->find('\'') !=std::string::npos)
                {
                    StrTokenize(*iBl, tBuf, '\'');
                }
                else
                {
                    StrTokenize(*iBl, tBuf);
                }
                if ((int)tBuf.size() > 1)
                {
                    for (unsigned i=1; i < tBuf.size(); i++ )
                    {
                        itsCurCryst->itsSpaceGroup->sgSymb["Hall"].push_back(tBuf[i]);
                    }
                }
                else
                {
                    std::cout << "Wrong Line format for the line contain '_symmetry_space_group_name_Hall' " 
                              << *iBl << std::endl;
                    errMsg.push_back("Wrong Line format for the line contain '_symmetry_space_group_name_Hall' ");
                    lErr = true;
                }
                
            }
            else if (iBl->find("_symmetry_space_group_name_H-M") !=std::string::npos)
            {
                std::vector<std::string> tBuf;
                if (iBl->find('\'') !=std::string::npos)
                {
                    StrTokenize(*iBl, tBuf, '\'');
                }
                else
                {
                    StrTokenize(*iBl, tBuf);
                }
                if ((int)tBuf.size() > 1)
                {
                    for (unsigned i=1; i < tBuf.size(); i++ )
                    {
                        itsCurCryst->itsSpaceGroup->sgSymb["xHM"].push_back(tBuf[i]);
                    }
                }
                else
                {
                    std::cout << "Wrong Line format for xHM symbol in line " << *iBl << std::endl;
                    errMsg.push_back("Wrong Line format for xHM symbol in line\n");
                    lErr = true;
              
                }
            }
            else if (iBl->find("_cell_length_a") !=std::string::npos)
            {
                std::vector<std::string> tBuf;
                StrTokenize(*iBl, tBuf);
                if ((int)tBuf.size() == 2)
                {
                    
                    if(tBuf[1].find("(") !=std::string::npos)
                    {
                        std::vector<ID> tBuf2;
                        StrTokenize(tBuf[1], tBuf2); 
                        itsCurCryst->itsCell->a = StrToReal(tBuf2[0]);
                    }
                    else
                    {
                        itsCurCryst->itsCell->a = StrToReal(tBuf[1]);
                    }
                }
                else
                {
                    std::cout << "No cell parameter a in line " << *iBl << std::endl;
                    errMsg.push_back("No cell parameter a in line \n");
                    lErr = true;
                }
            }
            else if (iBl->find("_cell_length_b") !=std::string::npos)
            {
                std::vector<std::string> tBuf;
                StrTokenize(*iBl, tBuf);
                if ((int)tBuf.size() == 2)
                {
                    if(tBuf[1].find("(") !=std::string::npos)
                    {
                        std::vector<ID> tBuf2;
                        StrTokenize(tBuf[1], tBuf2); 
                        itsCurCryst->itsCell->b = StrToReal(tBuf2[0]); 
                    }
                    else
                    {
                        itsCurCryst->itsCell->b = StrToReal(tBuf[1]);
                    }
                }
                else
                {
                    std::cout << "No space cell parameter b in line " << *iBl << std::endl;
                    errMsg.push_back("No cell parameter b in line \n");
                    lErr = true;
                }
            }
            else if (iBl->find("_cell_length_c") !=std::string::npos)
            {
                std::vector<std::string> tBuf;
                StrTokenize(*iBl, tBuf);
                if ((int)tBuf.size() == 2)
                {
                    if(tBuf[1].find("(") !=std::string::npos)
                    {
                        std::vector<ID> tBuf2;
                        StrTokenize(tBuf[1], tBuf2); 
                        itsCurCryst->itsCell->c = StrToReal(tBuf2[0]); 
                    }
                    else
                    {
                        itsCurCryst->itsCell->c = StrToReal(tBuf[1]);
                    }
                }
                else
                {errMsg.push_back("No cell parameter c in line \n");
                    lErr = true;
                    std::cout << "No cell parameter c in line " << *iBl << std::endl;
                    
                }
            }
            else if (iBl->find("_cell_angle_alpha") !=std::string::npos)
            {
                std::vector<std::string> tBuf;
                StrTokenize(*iBl, tBuf);
                if ((int)tBuf.size() == 2)
                {
                    
                    if(tBuf[1].find("(") !=std::string::npos)
                    {
                        std::vector<ID> tBuf2;
                        StrTokenize(tBuf[1], tBuf2, '('); 
                        itsCurCryst->itsCell->alpha = StrToReal(tBuf2[0]);
                    }
                    else
                    {
                        itsCurCryst->itsCell->alpha = StrToReal(tBuf[1]);
                    }
                }
                else
                {
                    std::cout << "No cell parameter alpha in line " << *iBl << std::endl;
                    errMsg.push_back("No cell parameter alpha in line \n");
                    lErr = false;
                }
            }
            else if (iBl->find("_cell_angle_beta") !=std::string::npos)
            {
                std::vector<std::string> tBuf;
                StrTokenize(*iBl, tBuf);
                if ((int)tBuf.size() == 2)
                {
                    
                    if(tBuf[1].find("(") !=std::string::npos)
                    {
                        std::vector<ID> tBuf2;
                        StrTokenize(tBuf[1], tBuf2, '('); 
                        itsCurCryst->itsCell->beta = StrToReal(tBuf2[0]);
                    }
                    else
                    {
                        itsCurCryst->itsCell->beta = StrToReal(tBuf[1]);
                    }
                }
                else
                {
                    std::cout << "No cell parameter beta in line " << *iBl << std::endl;
                    
                    errMsg.push_back("No cell parameter beta in line \n");
                    lErr = true;
                }
            }
            else if (iBl->find("_cell_angle_gamma") !=std::string::npos)
            {
                std::vector<std::string> tBuf;
                StrTokenize(*iBl, tBuf);
                if ((int)tBuf.size() == 2)
                {
                    
                    if(tBuf[1].find("(") !=std::string::npos)
                    {
                        std::vector<ID> tBuf2;
                        StrTokenize(tBuf[1], tBuf2, '('); 
                        itsCurCryst->itsCell->gamma = StrToReal(tBuf2[0]);
                    }
                    else
                    {
                        itsCurCryst->itsCell->gamma = StrToReal(tBuf[1]);
                    }
                }
                else
                {
                    std::cout << "No cell parameter gamma in line " << *iBl << std::endl;
                    errMsg.push_back("No cell parameter gamma in line \n");
                    lErr = true;
                }
            }
            else if (iBl->find("_cell_volume") !=std::string::npos)
            {
                std::vector<std::string> tBuf;
                StrTokenize(*iBl, tBuf);
                if ((int)tBuf.size() == 2)
                {
                    if(tBuf[1].find("(") !=std::string::npos)
                    {
                        std::vector<ID> tBuf2;
                        StrTokenize(tBuf[1], tBuf2); 
                        itsCurCryst->itsCell->vol = StrToReal(tBuf2[0]); 
                    }
                    else
                    {
                        itsCurCryst->itsCell->vol = StrToReal(tBuf[1]);
                    }
                }
                else
                {
                    std::cout << "No cell volume in line " << *iBl << std::endl;
                    errMsg.push_back("No cell volume in line \n");
                    lErr = true;
                }
            }
        }
        
        // Check the content obtained.
        std::cout << "Inside the crystal information container generated, there are, "  
                  << std::endl;
        std::cout << "1. Space Group related : " << std::endl;
        std::cout << "Space group number ";
        if (itsCurCryst->itsSpaceGroup->sgNum==-1)
        {
            std::cout << "is not presented in the input file " << std::endl;
        }
        else
        {
            std::cout << itsCurCryst->itsSpaceGroup->sgNum << std::endl;
        }
        
        std::cout << "Number of space group names " 
                  << itsCurCryst->itsSpaceGroup->sgSymb.size()
                  << std::endl;
      
                
        if (itsCurCryst->itsSpaceGroup->sgSymb.size() !=0)
        {
            std::cout << "Space group names : " << std::endl;
            for (std::map<ID, std::vector<ID> >::iterator iSN=itsCurCryst->itsSpaceGroup->sgSymb.begin();
                    iSN != itsCurCryst->itsSpaceGroup->sgSymb.end(); iSN++)
            {
                std::cout << iSN->first << " <--> ";
                for (std::vector<ID>::iterator iSV=iSN->second.begin();
                        iSV != iSN->second.end(); iSV++)
                {
                    std::cout << *iSV << "  " ;
                }
                std::cout << std::endl;
            }
        }
        
        if(itsCurCryst->itsCell)
        {
            std::cout << "2. Cell parameters : " << std::endl;
            std::cout << "Cell length a " << itsCurCryst->itsCell->a << std::endl
                      << "Cell length b " << itsCurCryst->itsCell->b << std::endl
                      << "Cell length c " << itsCurCryst->itsCell->c << std::endl
                      << "Cell angle alpha " << itsCurCryst->itsCell->alpha << std::endl
                      << "Cell angle beta "  << itsCurCryst->itsCell->beta << std::endl
                      << "Cell angle gamma " << itsCurCryst->itsCell->gamma << std::endl;
            
            if (itsCurCryst->itsCell->vol !=0.0)
            {
                std::cout << "Cell volume " << itsCurCryst->itsCell->vol << std::endl;
            }
            if (!itsCurCryst->itsCell->lattice.empty())
            {
                std::cout << "Crystal lattice " << itsCurCryst->itsCell->lattice << std::endl;
            }
        }
        
    }
    
    void GenCifFile::getCifCrystInfo(std::map<std::string,std::string>& tRowProps, 
                                     std::map<int,std::map<ID,std::vector<std::string> > >& tColProps)
    {
        if (itsCurCryst ==NullPoint)
        {
            itsCurCryst = new CrystInfo(); 
        }
        
        if (itsCurCryst->itsSpaceGroup == NullPoint)
        {
            itsCurCryst->itsSpaceGroup = new SpaceGroupMember();
        }
        
        if (itsCurCryst->itsCell==NullPoint)
        {
            itsCurCryst->itsCell  = new Cell();
        }
        
        
        
        if(itsCurCryst->itsResolution==NullPoint)
        {
            itsCurCryst->itsResolution = new Resolution();
        }
        
        
        // Space group related 
        
        
        if (tRowProps.find("_space_group_IT_number") !=tRowProps.end())
        {
            itsCurCryst->itsSpaceGroup->sgNum = StrToInt(tRowProps["_space_group_IT_number"]);
        }
        if (tRowProps.find("_symmetry_cell_setting") !=tRowProps.end())
        {
            itsCurCryst->itsCell->lattice =tRowProps["_symmetry_cell_setting"];
        }
        if (tRowProps.find("_symmetry_space_group_name_Hall") !=tRowProps.end())
        {
            itsCurCryst->itsSpaceGroup->sgSymb["Hall"].push_back(tRowProps["_symmetry_space_group_name_Hall"]);
        }
        if (tRowProps.find("_symmetry_space_group_name_H-M") !=tRowProps.end())
        {
            itsCurCryst->itsSpaceGroup->sgSymb["xHM"].push_back(tRowProps["_symmetry_space_group_name_H-M"]);
        }
        
        // Cell parameters related 
        if (tRowProps.find("_cell_length_a") !=tRowProps.end())
        {
            itsCurCryst->itsCell->a =StrToReal(tRowProps["_cell_length_a"]);
        }
        else
        {
            std::cout << "No cell length a in the input cif file " << std::endl;
            errMsg.push_back("No cell length a in the input cif file\n");
            lErr = true;
        }
        
        if (tRowProps.find("_cell_length_b") !=tRowProps.end())
        {
            itsCurCryst->itsCell->b =StrToReal(tRowProps["_cell_length_b"]);
        }
        else
        {
            std::cout << "No cell length b in the input cif file " << std::endl;
            errMsg.push_back("No cell length b in the input cif file\n");
            lErr = true;
        }
        
        if (tRowProps.find("_cell_length_c") !=tRowProps.end())
        {
            itsCurCryst->itsCell->c =StrToReal(tRowProps["_cell_length_c"]);
        }
        else
        {
            std::cout << "No cell length c in the input cif file " << std::endl;
            errMsg.push_back("No cell length c in the input cif file\n");
            lErr = true;
        }
        
        if (tRowProps.find("_cell_angle_alpha") !=tRowProps.end())
        {
            if(tRowProps["_cell_angle_alpha"].find("(") !=std::string::npos)
            {
                std::vector<ID> tBuf2;
                StrTokenize(tRowProps["_cell_angle_alpha"], tBuf2, '('); 
                itsCurCryst->itsCell->alpha = StrToReal(tBuf2[0]);
            }
            else
            {
                itsCurCryst->itsCell->alpha =StrToReal(tRowProps["_cell_angle_alpha"]);
            }
            
        }
        else
        {
            std::cout << "No cell angle alpha in the input cif file " << std::endl;
            errMsg.push_back("No cell angle alpha in the input cif file\n");
            lErr = true;
        }
        
        if (tRowProps.find("_cell_angle_beta") !=tRowProps.end())
        {
            if(tRowProps["_cell_angle_beta"].find("(") !=std::string::npos)
            {
                std::vector<ID> tBuf2;
                StrTokenize(tRowProps["_cell_angle_beta"], tBuf2, '('); 
                itsCurCryst->itsCell->beta = StrToReal(tBuf2[0]);
            }
            else
            {
                itsCurCryst->itsCell->beta =StrToReal(tRowProps["_cell_angle_beta"]);
            }
            
        }
        else
        {
            std::cout << "No cell angle beta in the input cif file " << std::endl;
            errMsg.push_back("No cell angle beta in the input cif file\n");
            lErr = true;
        }
        
        if (tRowProps.find("_cell_angle_gamma") !=tRowProps.end())
        {
            if(tRowProps["_cell_angle_gamma"].find("(") !=std::string::npos)
            {
                std::vector<ID> tBuf2;
                StrTokenize(tRowProps["_cell_angle_gamma"], tBuf2, '('); 
                itsCurCryst->itsCell->gamma = StrToReal(tBuf2[0]);
            }
            else
            {
                itsCurCryst->itsCell->gamma =StrToReal(tRowProps["_cell_angle_gamma"]);
            }
        }
        else
        {
            std::cout << "No cell angle gamma in the input cif file " << std::endl;
            errMsg.push_back("No cell angle gamma in the input cif file\n");
            lErr = true;
        }
        
        if (tRowProps.find("_cell_volume") !=tRowProps.end())
        {
            if(tRowProps["_cell_volume"].find("(") !=std::string::npos)
            {
                std::vector<ID> tBuf2;
                StrTokenize(tRowProps["_cell_volume"], tBuf2, '('); 
                itsCurCryst->itsCell->vol = StrToReal(tBuf2[0]);
            }
            else
            {
                itsCurCryst->itsCell->vol =StrToReal(tRowProps["_cell_volume"]);
            }
        }
        /*
        for (std::map<ID, ID>::iterator iM=tRowProps.begin();
                iM != tRowProps.end(); iM++)
        {
            std::cout << "Property  " << iM->first 
                      << " : Value " << iM->second << std::endl;
        }
       */
        
        // Resolution related properties 
        if (checkResol)
        {
            itsCurCryst->itsResolution->setResolLimit(RESOLTHRESHOLD_U);
            
            if (tRowProps.find("_diffrn_reflns_theta_max") ==tRowProps.end())
            {
                // No way to decide the resolution
                resolOK = false;
                errMsg.push_back("UNDEFINITED RESOLUTION: no theta_max to decide the                      resolution of data\n");
            }
            else if(tRowProps.find("_diffrn_radiation_wavelength")==tRowProps.end())
            {
                // No way to decide the resolution
                resolOK = false;
                errMsg.push_back("UNDEFINITED RESOLUTION: no X-ray wave length to decide the resolution of data\n");
            }
            else
            {
                if (tRowProps.find("_diffrn_radiation_wavelength") !=tRowProps.end())
                {
                    //std::cout << "_diffrn_radiation_wavelength : " << tRowProps["_diffrn_radiation_wavelength"] << std::endl;
                    itsCurCryst->itsResolution->wavLen = StrToReal(tRowProps["_diffrn_radiation_wavelength"]);
                }
        
                if (tRowProps.find("_diffrn_reflns_theta_max") !=tRowProps.end())
                {
                    itsCurCryst->itsResolution->thetaMax = StrToReal(tRowProps["_diffrn_reflns_theta_max"]);
                }
        
                if (itsCurCryst->itsResolution->wavLen >0.0 
                    && fabs(itsCurCryst->itsResolution->thetaMax) > 1.0 )
                {
                    itsCurCryst->itsResolution->setResol();
                    if (itsCurCryst->itsResolution->dMax  
                         > itsCurCryst->itsResolution->resolLimit)
                    {
                        resolOK = false;
                        ID aErr= "LOW RESOLUTION: Resolution " 
                             + RealToStr(itsCurCryst->itsResolution->dMax)
                             + " is lower than threshold " 
                             + RealToStr(itsCurCryst->itsResolution->resolLimit) 
                             + "\n";
                        errMsg.push_back(aErr);
                    }
                }
                else
                {
                    resolOK = false;
                    ID aErr= 
                    "UNDEFINITED RESOLUTION: Theta_max or x-ray wavelength is wrong\n";
                    errMsg.push_back(aErr);
                }
            }
        }
        
        // Check the content obtained.
        std::cout << "Inside the crystal information container generated, there are, "  
                  << std::endl;
        std::cout << "1. Space Group related : " << std::endl;
        
        if (itsCurCryst->itsSpaceGroup->sgNum !=-1)
        {
            std::cout << "Space group number " 
                      << itsCurCryst->itsSpaceGroup->sgNum << std::endl;
        }
                
        if (itsCurCryst->itsSpaceGroup->sgSymb.size() !=0)
        {
            
            std::cout << "Space group names : " << std::endl;
            for (std::map<ID, std::vector<ID> >::iterator iSN=itsCurCryst->itsSpaceGroup->sgSymb.begin();
                    iSN != itsCurCryst->itsSpaceGroup->sgSymb.end(); iSN++)
            {
                std::cout << iSN->first << " <--> ";
                for (std::vector<ID>::iterator iSV=iSN->second.begin();
                        iSV != iSN->second.end(); iSV++)
                {
                    std::cout << *iSV << "  " ;
                }
                std::cout << std::endl;
            }
        }

        if(itsCurCryst->itsCell)
        {
            std::cout << "2. Cell parameters : " << std::endl;
            std::cout << "Cell length a " << itsCurCryst->itsCell->a << std::endl
                      << "Cell length b " << itsCurCryst->itsCell->b << std::endl
                      << "Cell length c " << itsCurCryst->itsCell->c << std::endl
                      << "Cell angle alpha " << itsCurCryst->itsCell->alpha << std::endl
                      << "Cell angle beta "  << itsCurCryst->itsCell->beta << std::endl
                      << "Cell angle gamma " << itsCurCryst->itsCell->gamma << std::endl;
            
            if (itsCurCryst->itsCell->vol !=0.0)
            {
                std::cout << "Cell volume " << itsCurCryst->itsCell->vol << std::endl;
            }
            if (!itsCurCryst->itsCell->lattice.empty())
            {
                std::cout << "Crystal lattice " << itsCurCryst->itsCell->lattice << std::endl;
            }
        }
        
        if (itsCurCryst->itsResolution)
        {
            
            if (itsCurCryst->itsResolution->lSet)
            {
                std::cout << "Information about the resolution. " << std::endl;
                std::cout << "THe x-ray wavelength is " 
                          <<  itsCurCryst->itsResolution->wavLen << std::endl
                          << "The max theta is " 
                          << itsCurCryst->itsResolution->thetaMax << std::endl
                          << "The resolution is " 
                          << itsCurCryst->itsResolution->dMax << std::endl;   
            }
            else
            {
                std::cout << "No validated information on the resolution in the cif file"
                          << std::endl;
            }
        }
         
    }
    
    void GenCifFile::getCifSymOps(std::vector<std::vector<std::string> >::iterator iBs)
    {   
        if (itsCurCryst==NullPoint)
        {
            itsCurCryst = new CrystInfo();
        }
        if (itsCurCryst->itsSpaceGroup==NullPoint)
        {
            itsCurCryst->itsSpaceGroup = new SpaceGroupMember();
        }
        for (std::vector<std::string>::iterator iBl=iBs->begin();
                iBl != iBs->end(); iBl++)
        {
            if ((*iBl).find("_symmetry")==std::string::npos)
            {
                cleanChar(*iBl, '\'');
                std::vector<ID> tBuf;
                StrTokenize(*iBl, tBuf,',');
                StrToSymmOps(tBuf, itsCurCryst->itsSpaceGroup->sgOp[*iBl]);
            }
        }
        
        // Check
        std::cout << "Symmetry operators for the lattice: " << std::endl; 
        for (std::map<std::string, std::vector<std::vector<REAL> > >::iterator 
                iOp =  itsCurCryst->itsSpaceGroup->sgOp.begin();
                iOp != itsCurCryst->itsSpaceGroup->sgOp.end(); iOp++ )
        {
            std::cout << "For operator  "    << iOp->first 
                      << ", its matrix is: "  << std::endl;
            for (std::vector<std::vector<REAL> >::iterator iMatRow=iOp->second.begin();
                    iMatRow !=iOp->second.end(); iMatRow++)
            {
                for (std::vector<REAL>::iterator iMatEle=iMatRow->begin();
                        iMatEle !=iMatRow->end(); iMatEle++)
                {
                    std::cout << *iMatEle << "    ";
                }
                std::cout << std::endl;
            }
        }
    }
    
    void GenCifFile::getCifSymOps(std::map<ID,std::vector<std::string> >  & tOnePropGroup)
    {
        if (itsCurCryst==NullPoint)
        {
            itsCurCryst = new CrystInfo();
        }
        if (itsCurCryst->itsSpaceGroup==NullPoint)
        {
            itsCurCryst->itsSpaceGroup = new SpaceGroupMember();
        }
        
        int aPos=getKeyWordPos("_symmetry_equiv_pos_as_xyz", 
                               tOnePropGroup["lab"]);
        if(aPos==-1)
        {
            aPos=getKeyWordPos("_space_group_symop_operation_xyz",
                               tOnePropGroup["lab"]);
        }
        std::cout << "symm opt pos " << aPos << std::endl;
        
        if(aPos !=-1)
        {
            for (std::vector<std::string>::iterator iOps=tOnePropGroup["cont"].begin();
                    iOps !=tOnePropGroup["cont"].end(); iOps++)
            {
                std::vector<std::string> tBuf0, tBuf1, tBuf2;
                StrLower((*iOps));
                if ((*iOps).find('\'') !=std::string::npos)
                {
                    
                   // cleanChar(tOnePropGroup["cont"][aPos], '\'');
                   StrTokenize((*iOps), tBuf0, '\'');
                   for (std::vector<std::string>::iterator iT=tBuf0.begin();
                           iT !=tBuf0.end(); iT++)
                   {
                       if(!(*iT).empty())
                       {
                           tBuf1.push_back((*iT));
                       }
                   }
                }
                else
                {
                    StrTokenize((*iOps), tBuf1);
                }
                
                for (unsigned i=0; i < tBuf1.size(); i++)
                {
                    std::cout  << tBuf1[i] << std::endl;
                }
                
               
                
                StrTokenize(tBuf1[aPos], tBuf2,',');
                
                for (unsigned i=0; i < tBuf2.size(); i++)
                {
                    tBuf2[i] = TrimSpaces(tBuf2[i]);
                    //std::cout << "tBuf2 " << i << "  " << tBuf2[i] << std::endl;
                }
                
                StrToSymmOps(tBuf2, itsCurCryst->itsSpaceGroup->sgOp[(*iOps)]);
            }
        } 
        
        // Check
        std::cout << "Symmetry operators for the lattice: " << std::endl; 
        for (std::map<std::string, std::vector<std::vector<REAL> > >::iterator 
                iOp =  itsCurCryst->itsSpaceGroup->sgOp.begin();
                iOp != itsCurCryst->itsSpaceGroup->sgOp.end(); iOp++ )
        {
            std::cout << "For operator  "    << iOp->first 
                      << ", its matrix is: "  << std::endl;
            for (std::vector<std::vector<REAL> >::iterator iMatRow=iOp->second.begin();
                    iMatRow !=iOp->second.end(); iMatRow++)
            {
                for (std::vector<REAL>::iterator iMatEle=iMatRow->begin();
                        iMatEle !=iMatRow->end(); iMatEle++)
                {
                    std::cout << *iMatEle << "    ";
                }
                std::cout << std::endl;
            }
        }
        
    }
    
    void GenCifFile::getCifAtomInfo(std::vector<std::vector<std::string> >::iterator iBs)
    {
        
        for (std::vector<std::string>::iterator iBl=iBs->begin();
                iBl != iBs->end(); iBl++)
        {
            std::vector<std::string> tF;
            StrTokenize(TrimSpaces(*iBl), tF);
            
            // check which entries exist in an atom
            if ((int)tF.size() ==1 && 
                tF[0].find("_atom_site_")!=std::string::npos)
                
            {
                //std::cout << tF[0] << std::endl;
                if(tF[0].find("_atom_site_label") !=std::string::npos)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("id",curBlockLine));
                    std::cout << "column index for atom id " << curBlockLine << std::endl;
                }
                else if(tF[0].find("_atom_site_type_symbol") !=std::string::npos)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("chemType",curBlockLine));
                    std::cout << "column index for atom chemType " << curBlockLine << std::endl;
                }
                else if(tF[0].find("_atom_site_fract_x") !=std::string::npos)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("fract_x",curBlockLine) );
                    std::cout << "column index for atom fractional x " << curBlockLine << std::endl;
                }
                else if(tF[0].find("_atom_site_fract_y") !=std::string::npos)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("fract_y",curBlockLine) );
                    // std::cout << curBlockLine << std::endl;
                }
                else if(tF[0].find("_atom_site_fract_z") !=std::string::npos)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("fract_z",curBlockLine) );
                    //std::cout << curBlockLine << std::endl;
                }
                else if(tF[0].find("_atom_site_U_iso_or_equiv") !=std::string::npos)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("isoB",curBlockLine));
                    //std::cout << "BBBB " << curBlockLine << std::endl;
                }
                else if(tF[0].find("_atom_site_occupancy") !=std::string::npos)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("ocp",curBlockLine));
                    //std::cout << curBlockLine << std::endl;
                }
                else if(tF[0].find("_atom_site_symmetry_multiplicity") !=std::string::npos
                        || tF[0].find("_atom_site_site_symmetry_multiplicity") !=std::string::npos)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("symm_mul",curBlockLine));
                    std::cout << "site_symm_mul " << curBlockLine << std::endl;
                }
                
                curBlockLine++;
            }
            else if ((int)tF.size() >2 && tF[0].find("#") ==std::string::npos)
            {
                //std::cout << *iBl << std::endl;
                
                itsCurAtom = new AtomDict();
                bool getCoords = false;
                itsCurAtom->seriNum = itsCurAtomSeriNum;
                itsCurAtomSeriNum ++;
                if (hasProps["atom"].find("id") != hasProps["atom"].end())
                {
                    itsCurAtom->existProps["id"] =  hasProps["atom"]["resName"];
                    cleanSymbol(tF[hasProps["atom"]["id"]], "\"");
                    itsCurAtom->id = TrimSpaces(tF[hasProps["atom"]["id"]]);
                    // std::cout << "Its ID: " << itsCurAtom->id << std::endl;
                }
                if (hasProps["atom"].find("chemType") != hasProps["atom"].end())
                {
                    itsCurAtom->existProps["chemType"] =  hasProps["atom"]["chemType"];
                    std::string tS1(TrimSpaces(tF[hasProps["atom"]["chemType"]]));
                    std::string tS2;
                    for (unsigned i =0; i < tS1.size(); i++)
                    {
                        if( !std::isdigit(tS1[i]))
                        {
                            tS2 = tS2+tS1[i];
                        }
                        else
                        {
                            break;
                        }
                    }
                    // itsCurAtom->chemType = TrimSpaces(tF[hasProps["atom"]["chemType"]]);
                    itsCurAtom->chemType = tS2;;
                    // std::cout << "Its chemType : " << itsCurAtom->chemType << std::endl;
                }
                if (hasProps["atom"].find("isoB") != hasProps["atom"].end())
                {
                    itsCurAtom->existProps["isoB"] = hasProps["atom"]["isoB"];
                    std::string tS = TrimSpaces(tF[itsCurAtom->existProps["isoB"]]);
                    if (tS.find("(") ==std::string::npos)
                    {
                        itsCurAtom->isoB = StrToReal(tS);
                    }
                    else
                    {
                        std::vector<std::string> tSV;
                        StrTokenize(tS, tSV, '(');
                        if ((int)tSV.size() > 1)
                        {
                            itsCurAtom->isoB = StrToReal(tSV[0]);
                        }
                        else
                        {
                            std::cout << "string for iso_B is " << tS << std::endl;
                        }
                    }
                    //std::cout << "Here: Its U_iso : " << itsCurAtom->isoB << std::endl;
                }
                //if (hasProps["atom"].find("parCharge") != hasProps["atom"].end())
                //{
                //    itsCurAtom->existProps["parChange"] = hasProps["atom"]["parCharge"];
                //    itsCurAtom->parCharge = StrToReal(tF[itsCurAtom->existProps["parChange"]]);
                //    //std::cout << "Its partialCharge :" 
                //    //        << itsCurAtom->parCharge
                //    //        << std::endl;
                //}
                if (hasProps["atom"].find("ocp") != hasProps["atom"].end())
                {
                    itsCurAtom->existProps["ocp"] = hasProps["atom"]["ocp"];
                    itsCurAtom->ocp = StrToReal(TrimSpaces(tF[itsCurAtom->existProps["ocp"]]));
                    std::cout << "Its occupancy : " << itsCurAtom->ocp << std::endl;
                }
                if (hasProps["atom"].find("symm_mul") != hasProps["atom"].end())
                {
                    itsCurAtom->existProps["symm_mul"] = hasProps["atom"]["symm_mul"];
                    itsCurAtom->symmMult = StrToInt(TrimSpaces(tF[itsCurAtom->existProps["symm_mul"]]));
                    std::cout << "Its symm_mul : " << itsCurAtom->symmMult << std::endl;
                }
                
                if (hasProps["atom"].find("fract_x") != hasProps["atom"].end())
                {
                    itsCurAtom->existProps["fract_x"] = hasProps["atom"]["fract_x"];
                    std::string tSX = TrimSpaces(tF[itsCurAtom->existProps["fract_x"]]);
                     
                    if (tSX.find("(") ==std::string::npos)
                    {
                        // std::cout << StrToReal(tSX) << std::endl;
                        itsCurAtom->fracCoords[0]=StrToReal(tSX);
                    }
                    else
                    {   
                        std::vector<std::string> tXV;
                        StrTokenize(TrimSpaces(tSX), tXV, '(');
                        if ((int)tXV.size() >1)
                        {
                            itsCurAtom->fracCoords[0]=StrToReal(tXV[0]);
                        }
                        else
                        {
                            std::cout << "the string for coordinate x is " << tSX << std::endl;
                        }
                    }
                   
                    //std::cout << "Its (fractional) coord x : " << itsCurAtom->fracCoords[0] << std::endl;
                    
                }
                if (hasProps["atom"].find("fract_y") != hasProps["atom"].end())
                {
                    itsCurAtom->existProps["fract_y"] = hasProps["atom"]["fract_y"];
                    std::string tSY = TrimSpaces(tF[itsCurAtom->existProps["fract_y"]]);
                    if (tSY.find("(") ==std::string::npos)
                    {
                        itsCurAtom->fracCoords[1]=StrToReal(tSY);
                        
                    }
                    else
                    {
                        std::vector<std::string> tYV;
                        StrTokenize(TrimSpaces(tSY), tYV, '(');
                        if ((int)tYV.size() >1)
                        {
                            itsCurAtom->fracCoords[1]=StrToReal(tYV[0]);
                        }
                        else
                        {
                            std::cout << "the string for coordinate y is " 
                                      << tSY << std::endl;
                        }
                    }
                    std::cout << "Its (fractional) coord y : " 
                              << itsCurAtom->fracCoords[1] << std::endl;
                    
                }
                if (hasProps["atom"].find("fract_z") != hasProps["atom"].end())
                {
                    itsCurAtom->existProps["fract_z"] = hasProps["atom"]["fract_z"];
                    std::string tSZ = TrimSpaces(tF[itsCurAtom->existProps["fract_z"]]);
                    if (tSZ.find("(") ==std::string::npos)
                    {
                        itsCurAtom->fracCoords[2]=StrToReal(tSZ);
                    }
                    else
                    {
                        std::vector<std::string> tZV;
                        StrTokenize(TrimSpaces(tSZ), tZV, '(');
                        if ((int)tZV.size() >1)
                        {
                            itsCurAtom->fracCoords[2]=StrToReal(tZV[0]);
                        }
                        else
                        {
                            std::cout << "the string for coordinate z is " 
                                      << tSZ << std::endl;
                        }
                    }
                    //std::cout << "Its (fractional) coord z : " 
                    //          << itsCurAtom->fracCoords[2] << std::endl;
                    
                    getCoords = true;
                    
                }
               
                if (itsCurAtom->existProps.find("fract_x") 
                    !=itsCurAtom->existProps.end()
                    && getCoords)
                {
                    // TranslateIntoUnitCell(itsCurAtom->fracCoords);
                    FractToOrtho(itsCurAtom->fracCoords, itsCurAtom->coords, 
                                 itsCurCryst->itsCell->a, 
                                 itsCurCryst->itsCell->b, 
                                 itsCurCryst->itsCell->c, 
                                 itsCurCryst->itsCell->alpha,
                                 itsCurCryst->itsCell->beta, 
                                 itsCurCryst->itsCell->gamma);
                    hasCoords = true;
                }
                else if (itsCurAtom->existProps.find("x") !=itsCurAtom->existProps.end())
                {
                    OrthoToFract(itsCurAtom->coords, itsCurAtom->fracCoords, itsCurCryst->itsCell->a, 
                            itsCurCryst->itsCell->b, itsCurCryst->itsCell->c, itsCurCryst->itsCell->alpha,
                            itsCurCryst->itsCell->beta, itsCurCryst->itsCell->gamma);
                    // TranslateIntoUnitCell(itsCurAtom->fracCoords);
                    itsCurAtom->coords.clear();
                    FractToOrtho(itsCurAtom->fracCoords, itsCurAtom->coords, itsCurCryst->itsCell->a, 
                            itsCurCryst->itsCell->b, itsCurCryst->itsCell->c, itsCurCryst->itsCell->alpha,
                            itsCurCryst->itsCell->beta, itsCurCryst->itsCell->gamma);
                    hasCoords = true;
                }
                
                if (hasProps["atom"].find("chemType") == hasProps["atom"].end() 
                    && hasProps["atom"].find("id") != hasProps["atom"].end() )
                {
                    fromIdToChemType(itsCurAtom->id, itsCurAtom->chemType);
                }
               
                allAtoms.push_back(*itsCurAtom);
                delete itsCurAtom;
                itsCurAtom = NULL;
            }  
        }
        
    }
    
    void GenCifFile::getCifAtomInfo(std::map<ID,std::vector<std::string> >  & tOnePropGroup)
    {
        
        int pos1, pos2, pos3, pos4, pos5, posOcp, posCalc, posB;
        pos1 = getKeyWordPos("_atom_site_type_symbol", 
                               tOnePropGroup["lab"]);
        pos2 = getKeyWordPos("_atom_site_label", tOnePropGroup["lab"]);
        pos3 = getKeyWordPos("_atom_site_fract_x", 
                               tOnePropGroup["lab"]);
        pos4 = getKeyWordPos("_atom_site_fract_y", 
                               tOnePropGroup["lab"]);
        pos5 = getKeyWordPos("_atom_site_fract_z", 
                               tOnePropGroup["lab"]);
        
        posOcp = getKeyWordPos("_atom_site_occupancy", tOnePropGroup["lab"]);
        
        posCalc = getKeyWordPos("_atom_site_calc_flag", tOnePropGroup["lab"]);
        
        posB    = getKeyWordPos("_atom_site_U_iso_or_equiv", tOnePropGroup["lab"]);
        

        if ((pos1 !=-1 || pos2 != -1) && pos3 !=-1
             && pos4 !=-1 && pos5 !=-1)
        {
            // All required properties are there 
            for (std::vector<std::string>::iterator iAtm=tOnePropGroup["cont"].begin();
                    iAtm!=tOnePropGroup["cont"].end(); iAtm++)
            {
                std::vector<std::string> tBuf;
                StrTokenize((*iAtm), tBuf);
                if (tBuf.size()==tOnePropGroup["lab"].size())
                {
                    bool lCal = false;
                    if (posCalc !=-1)
                    {
                        if (tBuf[posCalc].find("dum") !=std::string::npos)
                        {
                            lCal = true;
                        }
                    }
                    if (!lCal)
                    {
                        getAtomInfoFromLine(tBuf, pos1, pos2, pos3, pos4, pos5, posOcp, posCalc, posB);
                    }
                }
                else
                {
                    std::cout << "Line: " << (*iAtm) << std::endl 
                             << " is not consistent with the following tags of this loop_ " 
                             << std::endl;
                    for (std::vector<std::string>::iterator iLab=tOnePropGroup["lab"].begin();
                    iLab!=tOnePropGroup["lab"].end(); iLab++)
                    {
                        std::cout << (*iLab) << std::endl;
                    }
                    errMsg.push_back("Line is not consistent with the following tags of this loop_\n");
                    lErr = true;
                  
                }
            }
            
            // Convert element symbols in cif from low cases to the proper ones 
            // (like 7037562.cif)
            
        }
        
    }
   
    void GenCifFile::getAtomInfoFromLine(std::vector<std::string> & tStrs,
                                         int tP1, int tP2, int tP3, int tP4, 
                                         int tP5, int tPOcp, int tPCalc, int tPB)
    {
        if (itsCurAtom !=NullPoint)
        {
            delete itsCurAtom;
            itsCurAtom = NullPoint;
        }
        itsCurAtom = new AtomDict();
        itsCurAtom->seriNum = itsCurAtomSeriNum;
        itsCurAtom->fromOrig = itsCurAtom->seriNum;
        itsCurAtom->isInPreCell = true;
        itsCurAtomSeriNum ++;

        if (tP1 !=-1)
        {
            ID tS1 = TrimSpaces(tStrs[tP1]);
            itsCurAtom->sMolType = tS1;
            itsCurAtom->chemType = tS1;
        }
        //std::cout << "Its chemType : " << itsCurAtom->chemType << std::endl;
        
        if (tP2 !=-1)
        {
            cleanSymbol(tStrs[tP2], "\"");
            itsCurAtom->id = TrimSpaces(tStrs[tP2]);
        }
        //std::cout << "Its ID: " << itsCurAtom->id << std::endl;
        
        std::string tSX = TrimSpaces(tStrs[tP3]);
                     
        if (tSX.find("(") ==std::string::npos)
        {
            // std::cout << StrToReal(tSX) << std::endl;
            itsCurAtom->fracCoords[0]=StrToReal(tSX);
        }
        else
        {  
            std::vector<std::string> tXV;
            StrTokenize(tSX, tXV, '(');
            if ((int)tXV.size() >1)
            {
                itsCurAtom->fracCoords[0]=StrToReal(tXV[0]);
                
            }
            else
            {
                std::cout << "the string for coordinate x is " << tSX << std::endl;
                errMsg.push_back("the string for coordinate x is wrong\n");
                lErr = true; 
            }
        }    
        //std::cout << "Its (fractional) coord x : " 
        //          << itsCurAtom->fracCoords[0] << std::endl;   

        std::string tSY = TrimSpaces(tStrs[tP4]);
                     
        if (tSY.find("(") ==std::string::npos)
        {
            itsCurAtom->fracCoords[1]=StrToReal(tSY);
        }
        else
        {  
            std::vector<std::string> tXV;
            StrTokenize(tSY, tXV, '(');
            if ((int)tXV.size() >1)
            {
                itsCurAtom->fracCoords[1]=StrToReal(tXV[0]);
                
            }
            else
            {
                std::cout << "the string for coordinate y is " << tSY << std::endl;
                errMsg.push_back("the string for coordinate y is wrong\n");
                lErr = true; 
              
            }
        }    
        //std::cout << "Its (fractional) coord y : " 
        //          << itsCurAtom->fracCoords[1] << std::endl;   
        
        std::string tSZ = TrimSpaces(tStrs[tP5]);
                     
        if (tSZ.find("(") ==std::string::npos)
        {
            itsCurAtom->fracCoords[2]=StrToReal(tSZ);
        }
        else
        {  
            std::vector<std::string> tXV;
            StrTokenize(tSZ, tXV, '(');
            if ((int)tXV.size() >1)
            {
                itsCurAtom->fracCoords[2]=StrToReal(tXV[0]);
                
            }
            else
            {
                std::cout << "the string for coordinate z is " << tSZ << std::endl;
                errMsg.push_back("the string for coordinate z is wrong\n");
                lErr = true;      
            }
        }    
        // std::cout << "Its (fractional) coord z : " 
        //          << itsCurAtom->fracCoords[2] << std::endl; 
        
        if (tPOcp > -1 && tPOcp < (int)tStrs.size())
        {
            std::string sOcp = TrimSpaces(tStrs[tPOcp]);
                     
            if (sOcp.find("(") ==std::string::npos)
            {
                if (sOcp.find("?") !=std::string::npos)
                {
                    itsCurAtom->ocp = 1.0;
                }
                else
                {
                    itsCurAtom->ocp=StrToReal(sOcp);
                }
            }
            else 
            {
                std::vector<std::string> tXV;
                StrTokenize(sOcp, tXV, '(');
                if ((int)tXV.size() >1)
                {
                    itsCurAtom->ocp=StrToReal(tXV[0]);
                }
                else
                {
                    std::cout << "Bug: The line related to Occupancy is-> "
                              << sOcp << std::endl;
                }
            }
            
            if (tPCalc > -1 && tPCalc < (int)tStrs.size())
            {
                std::string sCalc = TrimSpaces(tStrs[tPCalc]);
                if (sCalc.find("c") != std::string::npos)
                {
                    itsCurAtom->fromCalc = true;
                }
                else
                {
                    itsCurAtom->fromCalc = false;
                }
                
            }
            
        }    
        
        if (tPB > -1 && tPB < (int)tStrs.size())
        {

            itsCurAtom->isoB=StrToReal(cleanBrackets(tStrs[tPB]));  
        }
        
       
        FractToOrtho(itsCurAtom->fracCoords, itsCurAtom->coords, itsCurCryst->itsCell->a, 
                     itsCurCryst->itsCell->b, itsCurCryst->itsCell->c, itsCurCryst->itsCell->alpha,
                     itsCurCryst->itsCell->beta, itsCurCryst->itsCell->gamma);
        hasCoords = true;
       
        
        allAtoms.push_back(*itsCurAtom);
        delete itsCurAtom;
        itsCurAtom = NULL;
        
    }
    
    void GenCifFile::getCifAtomOxiInfo(std::map<ID,std::vector<std::string> >& tOnePropGroup)
    {
        int posType, posOxi;
        posType = getKeyWordPos("_atom_type_symbol", 
                               tOnePropGroup["lab"]);
        posOxi  = getKeyWordPos("_atom_type_oxidation_number", 
                                tOnePropGroup["lab"]);
        if (posType !=-1 && posOxi !=-1)
        {
            // All required properties are there 
            for (std::vector<std::string>::iterator iAtm=tOnePropGroup["cont"].begin();
                    iAtm!=tOnePropGroup["cont"].end(); iAtm++)
            {
                std::vector<std::string> tBuf;
                StrTokenize((*iAtm), tBuf);
                if (tBuf.size()==tOnePropGroup["lab"].size() 
                    && posType < tBuf.size() && posOxi < tBuf.size())
                {
                    getAtomOxiInfoFromLine(tBuf, posType, posOxi);
                }
                else
                {
                    std::cout << "Line: " << (*iAtm) << std::endl 
                             << " is not consistent with the following tags of this loop_ " 
                             << std::endl;
                    for (std::vector<std::string>::iterator iLab=tOnePropGroup["lab"].begin();
                    iLab!=tOnePropGroup["lab"].end(); iLab++)
                    {
                        std::cout << (*iLab) << std::endl;
                    }
                    errMsg.push_back("Line: is not consistent with the following tags of this loop_ "); 
                }
            }
        }
        
    }
    
    void GenCifFile::getAtomOxiInfoFromLine(std::vector<std::string>& tStrs, 
                                            int tPosType, int tPosOxi)
    {
        if (tStrs[tPosOxi].find("?") ==std::string::npos)
        {
            typeChargeMap[tStrs[tPosType]] = StrToReal(tStrs[tPosOxi]);
        }
    }
    
    void GenCifFile::setAtomOxiType()
    {
        if (allAtoms.size() !=0 && typeChargeMap.size() !=0)
        {
            if (allAtoms[0].sMolType != NullString)
            {
                for (std::vector<AtomDict>::iterator iAt= allAtoms.begin();
                        iAt != allAtoms.end(); iAt++)
                {
                    if (typeChargeMap.find(iAt->sMolType) != typeChargeMap.end())
                    {
                        iAt->formalCharge = typeChargeMap[iAt->sMolType];
                    }
                    else   
                    {
                        std::cout << "Check atom " << iAt->id 
                                  << " atom_site_type_symbol " << std::endl; 
                    }
                }
            }
        }
    }
    
    void GenCifFile::getHAtomIdxs()
    {
        allHAtomIdx.clear();
        
        for (std::vector<AtomDict>::iterator iAt=allAtoms.begin();
                iAt !=allAtoms.end(); iAt++)
        {
            if (iAt->chemType.compare("H") ==0)
            {
                allHAtomIdx.push_back(iAt->seriNum);
            }
        }
        
        if(allHAtomIdx.size() > 0)
        {
            hasH = true;
        }
    }
    
    
    
    void GenCifFile::reOrdErrMsg()
    {
        errMsg.clear();
        
        if (!notPowder)
        {
            errMsg.push_back("REJECTED STRUCTURE: Powder Diffraction ");
        }
        else if (!notNeuD)
        {
            errMsg.push_back("REJECTED STRUCTURE: Neutron Diffraction ");
        }
        else if (!resolOK)
        {
            if (checkResol)
            {
                std::string aMsg = "REJECTED STRUCTURE: Resolutions \n";
            
                
                if (itsCurCryst->itsResolution->dMax > 0.0)
                {
                    aMsg.append("LOW RESOLUTION: ");
                    aMsg.append(RealToStr(itsCurCryst->itsResolution->dMax)
                        + " is lower than threshold " 
                        + RealToStr(itsCurCryst->itsResolution->resolLimit) 
                        + "\n");  
                }
                else
                {
                    aMsg.append("(1) UNDEFINITED RESOLUTION: ");
                    aMsg.append("no theta_max to decide the resolution of data\n");
                    aMsg.append("or (2) UNDEFINITED RESOLUTION:  ");
                   aMsg.append("no X-ray wave length to decide the resolution of data\n");
                }
                errMsg.push_back(aMsg);
            }
        }
        else if (!RFactorOK)
        {
            std::string aMsg = "REJECTED STRUCTURE:  R factor related.\n";
            aMsg.append("(1) high R factors\n");
            aMsg.append("or (2) no R factors in the data\n");
            errMsg.push_back(aMsg);
        }
        else if (!colidOK)
        {
            errMsg.push_back("REJECTED STRUCTURE: Non-H atoms have less than 1.0 occp.\n");
        }
        else if (hasHeavyCalcAtoms)
        {
            errMsg.push_back("REJECTED STRUCTURE: Atoms of theoretically calculated found.\n");
        }
        else if (!symmOprOK)
        {
            errMsg.push_back("REJECTED STRUCTURE: No symmetry operators, can not generator molecules \n");
        }
        else if (lErr)
        {
            errMsg.push_back("REJECTED STRUCTURE: Reasons unknown\n");
        }
    }
    
    void GenCifFile::outAtomElems(Name tUserOutRoot)
    {
        Name initAtmElemsFName(tUserOutRoot);
        initAtmElemsFName.append("_init_atoms.list");
        std::map<std::string, std::vector<std::string> > initAtmElems;
        for (std::vector<AtomDict>::iterator iA = allAtoms.begin();
                            iA != allAtoms.end(); iA++)
        {
            initAtmElems[iA->chemType].push_back(iA->id);
        }
        
        if (initAtmElems.size() >0)
        {
            std::ofstream  initAtmElemsFile(initAtmElemsFName.c_str());
            if (initAtmElemsFile.is_open())
            {
                initAtmElemsFile.width(30);
                initAtmElemsFile << std::left << "ELEMENTS IN ATOM  : ";
                for (std::map<std::string, std::vector<std::string> >
                        ::iterator iE=initAtmElems.begin();
                        iE != initAtmElems.end(); iE++)
                {
                    initAtmElemsFile.width(8);
                    initAtmElemsFile  << std::left << iE->first;
                }
                initAtmElemsFile << std::endl;
                
                for (std::map<std::string, std::vector<std::string> >
                        ::iterator iE=initAtmElems.begin();
                        iE != initAtmElems.end(); iE++)
                {
                    initAtmElemsFile.width(30);
                    initAtmElemsFile << std::left << "Atoms with  "
                                     << iE->first  << "    :    ";
                    for (std::vector<std::string>::iterator 
                         iAId=iE->second.begin(); 
                         iAId != iE->second.end(); iAId++ )
                    {
                        initAtmElemsFile.width(8);
                        initAtmElemsFile << std::left << *iAId;
                    }
                    initAtmElemsFile << std::endl;   
                }
                
                initAtmElemsFile.close();
            }
        }
        
    }
    
    // ################################## another class for cif files of CCP4 
    // dictionary 
    
    
    
    DictCifFile::DictCifFile() : curBlockLine(ZeroInt),
            isPeptide(false),
            isDRna(false),
            hasConnect(false),
            hasCoords(false),
            hasH(false),
            itsCurAtomSeriNum(ZeroInt),
            itsCurAtom(NullPoint),
            itsCurBondSeriNum(ZeroInt),
            itsCurBond(NullPoint),
            itsCurAngleSeriNum(ZeroInt),
            itsCurAngle(NullPoint),
            itsCurTorsionSeriNum(ZeroInt),
            itsCurTorsion(NullPoint),
            itsCurChiralSeriNum(ZeroInt),
            itsCurChiral(NullPoint)
    {
    }
    
    DictCifFile::DictCifFile(Name tFname, 
                             std::ios_base::openmode tOpenMode) :
                             curBlockLine(ZeroInt),
                             isPeptide(false),
                             isDRna(false),
                             hasConnect(false),
                             hasCoords(false),
                             hasH(false),
                             hasCCP4Type(false),
                             itsCurAtomSeriNum(ZeroInt),
                             itsCurAtom(NullPoint),
                             itsCurBondSeriNum(ZeroInt),
                             itsCurBond(NullPoint),
                             itsCurAngleSeriNum(ZeroInt),
                             itsCurAngle(NullPoint),
                             itsCurTorsionSeriNum(ZeroInt),
                             itsCurTorsion(NullPoint),
                             itsCurChiralSeriNum(ZeroInt),
                             itsCurChiral(NullPoint)
    {
        
        if (tOpenMode == std::ios::in)
        {
            inFile.open(tFname.c_str(), tOpenMode);
            if (inFile.is_open())
            {
                itsCurAtom    = new AtomDict();
                itsCurBond    = new BondDict();
                itsCurAngle   = new AngleDict();
                itsCurTorsion = new TorsionDict();
                itsCurChiral  = new ChiralDict();
            
                setupSystem();
                 
            }
            else
            {
                std::cout << tFname << " can not be open for reading. Check the file "
                         << std::endl;
                exit(1);
            }
            
        }
        else
        {
            outFile.open(tFname.c_str(), tOpenMode);
        }
    }
    
    DictCifFile::DictCifFile(FileName                tFname,
                             std::ios::openmode      tOpenMode=std::ios::in ) :
                             curBlockLine(ZeroInt),
                             isPeptide(false),
                             isDRna(false),
                             hasConnect(false),
                             hasCoords(false),
                             hasH(false),
                             hasCCP4Type(false),
                             itsCurAtomSeriNum(ZeroInt),
                             itsCurAtom(NullPoint),
                             itsCurBondSeriNum(ZeroInt),
                             itsCurBond(NullPoint),
                             itsCurAngleSeriNum(ZeroInt),
                             itsCurAngle(NullPoint),
                             itsCurTorsionSeriNum(ZeroInt),
                             itsCurTorsion(NullPoint),
                             itsCurChiralSeriNum(ZeroInt),
                             itsCurChiral(NullPoint)
    {   
        propComp.id       = NullString;
        propComp.code     = NullString;
        propComp.name     = NullString;
        propComp.group    = NullString;
        propComp.numAtoms = 0;
        propComp.numH     = 0;
        
        if (tOpenMode == std::ios::in)
        {
            inFile.open(tFname, tOpenMode);
            if (inFile.is_open())
            {
                itsCurAtom    = new AtomDict();
                itsCurBond    = new BondDict();
                itsCurAngle   = new AngleDict();
                itsCurTorsion = new TorsionDict();
                itsCurChiral  = new ChiralDict();
                
                setupSystem();
            }
            else
            {
                std::cout << tFname << " can not be open for reading. Check the file ! "
                         << std::endl;
                exit(1);
            }
            
        }
        else
        {
            outFile.open(tFname, tOpenMode);
        }        
    }
    
    
    DictCifFile::DictCifFile(FileName                tCifName,
                             FileName                tPdbName) :
                             curBlockLine(ZeroInt),
                             isPeptide(false),
                             isDRna(false),
                             hasConnect(false),
                             hasCoords(false),
                             hasH(false),
                             hasCCP4Type(false),
                             itsCurAtomSeriNum(ZeroInt),
                             itsCurAtom(NullPoint),
                             itsCurBondSeriNum(ZeroInt),
                             itsCurBond(NullPoint),
                             itsCurAngleSeriNum(ZeroInt),
                             itsCurAngle(NullPoint),
                             itsCurTorsionSeriNum(ZeroInt),
                             itsCurTorsion(NullPoint),
                             itsCurChiralSeriNum(ZeroInt),
                             itsCurChiral(NullPoint)
    {   
        std::ifstream aInCif(tCifName);
        
        if (aInCif.is_open())
        {
            
            itsCurAtom    = new AtomDict();
            setupSystem3Secs(aInCif);  
            aInCif.close();
          
            
            DictPDBFile aInPdb(tPdbName, std::ios::in);
            transCoordsPdbToCif(aInPdb);
            
        }
        else
        {
            std::cout << tCifName << " can not be open for reading. Check the file ! "
                      << std::endl;
            exit(1);
        }
    }
    
    DictCifFile::~DictCifFile()
    {
        if(inFile.is_open())
        {
            inFile.close();
        }
        
        if(outFile.is_open())
        {
            outFile.close();
        }
        if (itsCurAtom)
        {
            delete itsCurAtom; 
            itsCurAtom = NULL;
        }
        if(itsCurBond)
        {
            delete itsCurBond;
            itsCurBond = NULL;
        }
        if(itsCurAngle)
        {
            delete itsCurAngle;
            itsCurAngle = NULL;
        }
        if(itsCurTorsion)
        {
            delete itsCurTorsion;
            itsCurTorsion = NULL;
        }
        if(itsCurChiral)
        {
            delete itsCurChiral;
            itsCurChiral = NULL;
        }
    }
    
    
    
    void DictCifFile::initCifKeys()
    {
        // could be read from a file currently set up  sections of data description 
        // atoms, bonds
        
        // data descreption section 
        allCifKeys["dataDesc"].push_back("id");
        allCifKeys["dataDesc"].push_back("name");
        allCifKeys["dataDesc"].push_back("type");
        allCifKeys["dataDesc"].push_back("group");
        allCifKeys["dataDesc"].push_back("pdbx_type");
        allCifKeys["dataDesc"].push_back("formula");
        allCifKeys["dataDesc"].push_back("mon_nstd_parent_comp_id");
        allCifKeys["dataDesc"].push_back("pdbx_synonyms");
        allCifKeys["dataDesc"].push_back("pdbx_formal_charge");
        allCifKeys["dataDesc"].push_back("pdbx_initial_date");
        allCifKeys["dataDesc"].push_back("pdbx_modified_date");
        allCifKeys["dataDesc"].push_back("pdbx_ambiguous_flag");
        allCifKeys["dataDesc"].push_back("pdbx_release_status");
        allCifKeys["dataDesc"].push_back("pdbx_replaced_by");
        allCifKeys["dataDesc"].push_back("pdbx_replaces");
        allCifKeys["dataDesc"].push_back("formula_weight");
        allCifKeys["dataDesc"].push_back("one_letter_code");
        allCifKeys["dataDesc"].push_back("three_letter_code");
        allCifKeys["dataDesc"].push_back("pdbx_model_coordinates_details");
        allCifKeys["dataDesc"].push_back("pdbx_model_coordinates_missing_flag");
        allCifKeys["dataDesc"].push_back("pdbx_ideal_coordinates_details");
        allCifKeys["dataDesc"].push_back("pdbx_ideal_coordinates_missing_flag");
        allCifKeys["dataDesc"].push_back("pdbx_model_coordinates_db_code");
        allCifKeys["dataDesc"].push_back("pdbx_subcomponent_list");
        allCifKeys["dataDesc"].push_back("pdbx_processing_site");
        
        // atom section, check the doc for cif format 
        // hopefully those are enough for dictionary and pdb source
        allCifKeys["atom"].push_back("comp_id");
        allCifKeys["atom"].push_back("atom_id");
        allCifKeys["atom"].push_back("alt_atom_id");
        allCifKeys["atom"].push_back("type_symbol");
        allCifKeys["atom"].push_back("type_energy");
        allCifKeys["atom"].push_back("charge");
        allCifKeys["atom"].push_back("partial_charge");
        allCifKeys["atom"].push_back("pdbx_align");
        allCifKeys["atom"].push_back("pdbx_aromatic_flag");
        allCifKeys["atom"].push_back("pdbx_leaving_atom_flag");
        allCifKeys["atom"].push_back("pdbx_stereo_config");
        allCifKeys["atom"].push_back("x");
        allCifKeys["atom"].push_back("y");
        allCifKeys["atom"].push_back("z");
        allCifKeys["atom"].push_back("model_Cartn_x");
        allCifKeys["atom"].push_back("model_Cartn_y");
        allCifKeys["atom"].push_back("model_Cartn_z");
        allCifKeys["atom"].push_back("pdbx_model_Cartn_x_ideal");
        allCifKeys["atom"].push_back("pdbx_model_Cartn_y_ideal");
        allCifKeys["atom"].push_back("pdbx_model_Cartn_z_ideal");
        allCifKeys["atom"].push_back("pdbx_component_atom_id");
        allCifKeys["atom"].push_back("pdbx_component_comp_id");
        allCifKeys["atom"].push_back("pdbx_ordinal");
        
        // bond section
        allCifKeys["bond"].push_back("comp_id");
        allCifKeys["bond"].push_back("atom_id_1");
        allCifKeys["bond"].push_back("atom_id_2");
        allCifKeys["bond"].push_back("value_order");
        allCifKeys["bond"].push_back("pdbx_aromatic_flag");
        allCifKeys["bond"].push_back("pdbx_stereo_config");
        allCifKeys["bond"].push_back("pdbx_ordinal");
        allCifKeys["bond"].push_back("type");
        allCifKeys["bond"].push_back("value_dist");
        allCifKeys["bond"].push_back("value_dist_esd");
        
    }
    
    void DictCifFile::setupSystem()
    {
        
        initCifKeys();
        
        if (inFile.is_open() )
        { 
            bool tOK       = true;
            
            std::map<std::string, bool>    lBloc;
            
            lBloc["head"]     = false;
            lBloc["compList"] = false;
            lBloc["compData"] = false;
            lBloc["dataDesc"] = false;
            
            lBloc["atom"]     = false;
            lBloc["bond"]     = false;
            lBloc["angle"]    = false;
            lBloc["torsion"]  = false;
            lBloc["chiral"]   = false;
            lBloc["metal"]    = false;
            
            
            std::string tRecord="";
            
            while(!inFile.eof() && tOK)
            {
                std::getline(inFile, tRecord);
                tRecord = TrimSpaces(tRecord);
                
                std::vector<std::string> tBuf;
                std::vector<std::string> tBuf_t;
                if (tRecord.find('\"') !=std::string::npos)
                {
                    // std::cout << tRecord << std::endl;
                    
                    StrTokenize(tRecord, tBuf_t, '\"');
                    
                    for (int i=0; i < (int)tBuf_t.size(); i++)
                    {
                        if ((i+1)%2 !=0)
                        {
                            std::vector<std::string> tBuf_t2;
                            StrTokenize(tBuf_t[i], tBuf_t2);
                            for (int j=0; j< (int)tBuf_t2.size(); j++)
                            {
                                tBuf.push_back(tBuf_t2[j]);
                            }
                        }
                        else
                        {
                            // entries within the reversed comma 
                            tBuf.push_back(tBuf_t[i]);
                        }
                    }
                }
                else
                {
                    StrTokenize(tRecord, tBuf);
                }
                
                
                if ((int)tBuf.size() !=0)
                {
                    checkBloc(lBloc, tBuf);
                    // std::cout << tRecord << std::endl;
                    if(lBloc["head"])
                    {
                        //std::cout << "get header info" << std::endl;
                        getHeadInfo(tBuf);
                    }
                    else if(lBloc["compList"])
                    {
                        // std::cout << "get compound info" << std::endl;
                        getChemInfo(tRecord);
                    }
                    else if (lBloc["dataDesc"])
                    {
                        getDataDescription(tBuf);
                    }
                    else if(lBloc["atom"])
                    {   
                        //std::cout << "get atom info" << std::endl;
                        getAtomInfo(tBuf);
                    }
                    else if(lBloc["bond"])
                    {
                        //std::cout << "get bond info" << std::endl;
                        //std::cout << tRecord << std::endl;
                        getBondAndConnect(tBuf);
                    }
                    else if(lBloc["metal"])
                    {
                        //std::cout << "get metal info" << std::endl;
                        //std::cout << tRecord << std::endl;
                        getMetalCNGeo(tBuf);
                    }
                    //else if(lBloc["angle"])
                    //{
                        // std::cout << "get angle info " << std::endl;
                        
                        // getAngleInfo(tBuf);
                    //}
                    //else if(lBloc["torsion"])
                    //{
                    //    getTorsionInfo(tBuf);
                    //}
                    else if (lBloc["chiral"])
                    {
                        // std::cout << tRecord << std::endl;
                        getChiralInfo(tBuf);
                    }
                }
            }
            
            inFile.close();
            
            
            /*
            //std::cout << "CCP4 type " << hasCCP4Type << std::endl;
            
            
            std::cout << "The following is the property of the system in the input cif: "
                      << std::endl;
            
            std::cout << "The system ID " << propComp.id << std::endl;
            std::cout << "The system three-letter code " << propComp.code << std::endl;
            std::cout << "The system name  "  << propComp.name  << std::endl;
            std::cout << "The system group " << propComp.group << std::endl;
            std::cout << "Number of atoms in the system  " << propComp.numAtoms << std::endl;
            std::cout << "Number of H atoms in the system " << propComp.numH << std::endl;
            std::cout << " Is it a peptide ? ";
            if (isPeptide)
            {
                std::cout << " Yes " << std::endl;
            }
            else
            {
                std::cout << " No " << std::endl;   
            }
            
            */
            // std::cout << "The system description level " << propComp.level << std::endl;
            
            // Set the bonding properties for all atoms based
            // on their connections
            // If no bonds defined in the cif file. Stop the program at the moment.
            if ((int)allBonds.size()==0)
            {
                std::cout << "There is no bond defined in the cif file. Program stopped"
                        << std::endl;
                exit(1);    
            }
            //std::cout << "Number of atoms " << allAtoms.size() << std::endl;
            //std::cout << "Number of bonds " << allBonds.size() << std::endl;
                      
            
            setHydroAtomConnect();
            // addMissHydroAtoms();
            
            setAtomsBondingAndChiralCenter(allAtoms);
            
            // setAllAngles();
            
            setAtomsCChemType();
            
            
            
            setAtomsMetalType();
            
            setAtomsVDWRadius();
            
            setAtomsPartialCharges();
           
            
            ringDetecting();
            
            setAtomFormTypes(allAtoms);
            
            
            if (!hasCCP4Type)
            {
                setAtomsCCP4Type();
            }
            
            
            
            for (std::vector<AtomDict>::iterator iA = allAtoms.begin();
                    iA != allAtoms.end(); iA++)
            {
                std::cout << "\nAtom " << iA->seriNum << " : " << std::endl
                        << "Its ID : " << iA->id << std::endl
                        << "Its Chemical Type : " << iA->chemType << std::endl
                        << "Its CCP4 chemical type " << iA->ccp4Type << std::endl
                        << "Its bonding index : "   << iA->bondingIdx << std::endl
                        << "Its residue Name: " << iA->resName<< std::endl
                        << "Its formal charge " << iA->parCharge << std::endl;
                std::cout << "Its connected atoms are : " << std::endl;
                for (std::vector<int>::iterator iSer= iA->connAtoms.begin();
                        iSer != iA->connAtoms.end(); iSer++)
                {
                    std::cout << allAtoms[*iSer].id << std::endl;
                }  
                
                if (hasCoords)
                {
                    std::cout << "The coords for the atom are : " << std::endl;
                    
                    for (std::vector<REAL>::iterator iX= iA->coords.begin();
                        iX !=iA->coords.end(); iX++)
                    {
                        std::cout << *iX << std::endl;
                    }
                }
                
            }
            
            
            for (std::vector<BondDict>::iterator iBo = allBonds.begin();
                    iBo != allBonds.end(); iBo++)
            {
                std::cout << "For Bond : " << iBo->seriNum << std::endl;
                std::cout << "It is in Residue: " << iBo->resName << std::endl;
                std::cout << "Its component atom1 " << iBo->atoms[0] << std::endl;
                std::cout << "Its component atom2 " << iBo->atoms[1] << std::endl;
                //std::cout << "Its length : "    << iBo->length << std::endl;
                //std::cout << "Its sigLength : " << iBo->sigLength << std::endl;
                std::cout << "its order : " << iBo->order << std::endl;
            } 
            
            
        }
        
    }
    
    void DictCifFile::setupSystem3Secs(std::ifstream & tInCif)
    {
        bool tOK       = true;
            
        std::map<std::string, bool>    lBloc;
            
        lBloc["head"]     = false;
        lBloc["atom"]     = false;
        lBloc["others"]   = false;
            
        std::string tRecord="";
            
        while(!inFile.eof() && tOK)
        {
            std::getline(tInCif, tRecord);
            tRecord = TrimSpaces(tRecord);
                
            std::vector<std::string> tBuf;
            std::vector<std::string> tBuf_t;
            if (tRecord.find('\"') !=std::string::npos)
            {
                // std::cout << tRecord << std::endl;
                StrTokenize(tRecord, tBuf_t, '\"');
                    
                for (int i=0; i < (int)tBuf_t.size(); i++)
                {
                    if ((i+1)%2 !=0)
                    {
                        std::vector<std::string> tBuf_t2;
                        StrTokenize(tBuf_t[i], tBuf_t2);
                        for (int j=0; j< (int)tBuf_t2.size(); j++)
                        {
                            tBuf.push_back(tBuf_t2[j]);
                        }
                    }
                    else
                    {
                        // entries within the reversed comma 
                        tBuf.push_back(tBuf_t[i]);
                    }
                }
            }
            else
            {
                StrTokenize(tRecord, tBuf);
            }
            
            if ((int)tBuf.size() !=0)
            {
                checkBloc(lBloc, tBuf);
                // std::cout << tRecord << std::endl;
                if(lBloc["head"])
                {
                    allUnchangedBlocks["head"].push_back(tRecord);
                }
                else if(lBloc["atom"])
                {   
                    getAtomInfo(tBuf);
                }
                else 
                {
                    allUnchangedBlocks["others"].push_back(tRecord);  
                }
            }
        }
    }
        
    void DictCifFile::checkBloc(std::map<std::string, bool> & tL, 
                       std::vector<std::string> tF)
    {
        if((int)tF[0].size() >=1 and (int)tF[0].find("#") !=0)
        {
            if ((int)tF.size() == 1)
            {
                // std::cout << tF[0] << std::endl;
                
                if (TrimSpaces(tF[0]).find("loop_") !=std::string::npos)
                {
                    setFlags(tL, "loop");
                }
                else if (TrimSpaces(tF[0]).compare("global_") ==0)
                {
                    //std::cout<< "Setup flag head " <<std::endl;
                    setFlags(tL, "head");
                }
                else if (! tL["compList"] && 
                         TrimSpaces(tF[0]).find("_chem_comp.") 
                         !=std::string::npos)
                {
                    //std::cout << "setup compound info " << std::endl;
                    curBlockLine = 0;
                    setFlags(tL, "compList");
                    
                }
                else if ( ! tL["atom"] && 
                         TrimSpaces(tF[0]).find("_chem_comp_atom")
                         !=std::string::npos)
                {
                    // std::cout<< "Setup flag atom " <<std::endl;
                    curBlockLine = 0;
                    setFlags(tL, "atom");
                }
                else if ( ! tL["bond"] && 
                         TrimSpaces(tF[0]).find("_chem_comp_bond")
                         !=std::string::npos)
                {
                    // std::cout << "setup bond flag " << std::endl;
                    curBlockLine = 0;
                    setFlags(tL, "bond");
                }
                else if (not tL["angle"] && 
                         TrimSpaces(tF[0]).find("_chem_comp_angle")
                         !=std::string::npos)
                {
                    //std::cout << "setup bond angle flag " << std::endl;
                    curBlockLine = 0;
                    setFlags(tL, "angle");
                }
                else if (not tL["torsion"] && 
                         TrimSpaces(tF[0]).find("_chem_comp_tor")
                         !=std::string::npos)
                {
                    curBlockLine = 0;
                    setFlags(tL, "torsion");
                }
                else if (not tL["chiral"] && 
                         TrimSpaces(tF[0]).find("_chem_comp_chir")
                         !=std::string::npos)
                {
                    curBlockLine = 0;
                    setFlags(tL, "chiral");
                }
                else if (not tL["metal"] && TrimSpaces(tF[0]).find("_chem_comp_metal") 
                        !=std::string::npos)
                {
                    std::cout << "setup metal flag " << std::endl;
                    curBlockLine = 0;
                    setFlags(tL, "metal");
                }       
            }
            // the following is from pdb ligand cif 
            else if ((int)tF.size() == 2)
            {
                if (TrimSpaces(tF[0]).find("_chem_comp.") 
                         !=std::string::npos)
                {
                    curBlockLine = 0;
                    setFlags(tL, "dataDesc");
                }
            }
        }
    }
    
    void DictCifFile::setFlags(std::map<std::string, bool> & tL,
                               std::string tS = NullString)
    {
        for (std::map<std::string, bool>::iterator aT=tL.begin();
                aT !=tL.end(); aT++)
        {
            if (aT->first.compare(tS) == 0)
            {
                aT->second = true;
            }
            else
            {
                aT->second = false;
            }
        }
    }
    
    void DictCifFile::getHeadInfo(std::vector<std::string> tF)
    {
        if ((int)tF.size() ==2)
        {
           // std::cout << tF[0] << "  " << tF[1] << std::endl;
           
           if(tF[0].find(".name") !=std::string::npos)
           {
               dictCifHead.libName = TrimSpaces(tF[1]);
           }
           else if (tF[0].find(".version") !=std::string::npos)
           {
               dictCifHead.libVersion = TrimSpaces(tF[1]);
           }
           else if (tF[0].find(".type") !=std::string::npos)
           {
               dictCifHead.monType = TrimSpaces(tF[1]);
           }
           else if (tF[0].find(".group") !=std::string::npos)
           {
               dictCifHead.group = TrimSpaces(tF[1]);
           }
           else if (tF[0].find(".update") !=std::string::npos)
           {
               dictCifHead.libUpdate = TrimSpaces(tF[1]);
           }
        }
    }
    
    void DictCifFile::getChemInfo(std::string tRe)
    {
        if ((int)tRe.size() > 0)
        {
            std::vector<std::string> tF;
            StrTokenize(tRe, tF);
            
            if ((int)tF.size() ==1 && 
                    tF[0].find("_chem_comp.")!=std::string::npos)
            {
                std::vector<std::string> tF1;
                tF1 = StrTokenize(TrimSpaces(tF[0]), '.');
                // std::cout <<  tF[0] << std::endl;
                
                if((int)tF1.size() ==2)
                {
                    if(tF1[1].find("id") !=std::string::npos)
                    {
                        hasProps["compoundInfo"].insert(std::pair<std::string, int>("id",curBlockLine) );
                        curBlockLine++;
                    }
                    else if(tF1[1].find("three_letter_code") !=std::string::npos)
                    {
                        hasProps["compoundInfo"].insert(std::pair<std::string, int>("three_letter_code",curBlockLine));
                        curBlockLine++;
                    }
                    else if(tF1[1].find("name") !=std::string::npos)
                    {
                        hasProps["compoundInfo"].insert(std::pair<std::string, int>("name",curBlockLine) );
                        curBlockLine++;
                    }
                    else if(tF1[1].find("type") !=std::string::npos)
                    {
                        hasProps["compoundInfo"].insert(std::pair<std::string, int>("type",curBlockLine) );
                        curBlockLine++;
                    }
                    else if(tF1[1].find("group") !=std::string::npos)
                    {
                        hasProps["compoundInfo"].insert(std::pair<std::string, int>("group",curBlockLine) );
                        curBlockLine++;
                    }
                    else if(tF1[1].find("number_atoms_all") !=std::string::npos)
                    {
                        hasProps["compoundInfo"].insert(std::pair<std::string, int>("number_atoms_all",curBlockLine));
                        curBlockLine++;
                    }
                    else if(tF1[1].find("number_atoms_nh") !=std::string::npos)
                    {
                        hasProps["compoundInfo"].insert(std::pair<std::string, int>("number_atoms_nh",curBlockLine));
                        curBlockLine++;
                    }
                    else if(tF1[1].find("desc_level") !=std::string::npos)
                    {
                        hasProps["compoundInfo"].insert(std::pair<std::string, int>("desc_level",curBlockLine));
                        curBlockLine++;
                    }
                    
                }
            }
            
            
            
            tF.clear();
            
            
            tF = StrTokenize(tRe, '\'');       
            
            if ((int)tF.size() == 3)
            {
                std::vector<std::string> tF1;
                std::vector<std::string> tF2;
             
                StrTokenize(TrimSpaces(tF[0]), tF1);
                StrTokenize(TrimSpaces(tF[2]), tF2);
                
                std::vector<std::string> tF4;
                for(std::vector<std::string>::iterator iS=tF1.begin();
                        iS != tF1.end(); iS++)
                {
                    tF4.push_back(*iS);
                }
                tF4.push_back("\'"+ tF[1] + "\'");
                for(std::vector<std::string>::iterator iS=tF2.begin();
                        iS != tF2.end(); iS++)
                {
                    tF4.push_back(*iS);
                }
                int tNumEntry = (int)hasProps["compoundInfo"].size();
                //std::cout << "n-entry-defined " << tNumEntry << std::endl;
                //std::cout << "number of value entries : " << (int)tF4.size() << std::endl;
                if ((int)tF4.size() ==tNumEntry )
                {
                    
                    if (hasProps["compoundInfo"].find("id") != hasProps["compoundInfo"].end())
                    {
                        propComp.id   = TrimSpaces(tF4[hasProps["compoundInfo"]["id"]]);
                        //std::cout << "Compound ID " << propComp.id  << std::endl;
                    }
                    if (hasProps["compoundInfo"].find("three_letter_code") != hasProps["compoundInfo"].end())
                    {
                        propComp.code = TrimSpaces(tF4[hasProps["compoundInfo"]["three_letter_code"]]);
                        //std::cout << "Compound Code " << propComp.code  << std::endl;
                    }
                    if (hasProps["compoundInfo"].find("name") != hasProps["compoundInfo"].end())
                    {
                        propComp.name = TrimSpaces(tF4[hasProps["compoundInfo"]["name"]]);
                        //std::cout << "Compound name " << propComp.name  << std::endl;
                    }
                    if (hasProps["compoundInfo"].find("type") != hasProps["compoundInfo"].end())
                    {
                        propComp.group = TrimSpaces(tF4[hasProps["compoundInfo"]["type"]]);
                        //std::cout << "Compound type " << propComp.group  << std::endl;
                        StrUpper(propComp.group);
                        if (propComp.group.find("PEPTIDE") !=std::string::npos)
                        {
                            isPeptide = true;
                        }
                        else if (propComp.group.find("DNA") !=std::string::npos
                                 || propComp.group.find("RNA") !=std::string::npos)
                        {
                            isDRna = true;
                        }
                    }
                    if (hasProps["compoundInfo"].find("group") != hasProps["compoundInfo"].end())
                    {
                        propComp.group = TrimSpaces(tF4[hasProps["compoundInfo"]["group"]]);
                        //std::cout << "Compound group " << propComp.group  << std::endl;
                        StrUpper(propComp.group);
                        if (propComp.group.find("PEPTIDE") !=std::string::npos)
                        {
                            isPeptide = true;
                        }
                        else if (propComp.group.find("DNA") !=std::string::npos
                                 || propComp.group.find("RNA") !=std::string::npos)
                        {
                            isDRna = true;
                        }
                    }
                    if (hasProps["compoundInfo"].find("number_atoms_all") != hasProps["compoundInfo"].end())
                    {
                        std::string tST1 = tF4[hasProps["compoundInfo"]["number_atoms_all"]];
                        if (isInt(tST1))
                        {
                            propComp.numAtoms = StrToInt(tST1);
                            //std::cout << "Number of atoms  " << propComp.numAtoms  << std::endl;
                        }
                    }
                    if (hasProps["compoundInfo"].find("number_atoms_nh") != hasProps["compoundInfo"].end())
                    {
                        std::string tST2 = tF4[hasProps["compoundInfo"]["number_atoms_nh"]];
                        if (isInt(tST2))
                        {
                            propComp.numH = StrToInt(tST2);
                            //std::cout << "number of H atoms " << propComp.numH  << std::endl;
                        }
                    }
                    if (hasProps["compoundInfo"].find("desc_level") != hasProps["compoundInfo"].end())
                    {
                        propComp.level = TrimSpaces(tF4[hasProps["compoundInfo"]["desc_level"]]);
                        // std::cout << "Description level " << propComp.level << std::endl;
                    }
                }
            }
        }
    }
    
    void DictCifFile::getDataDescription(std::vector<std::string> tF)
    {
        std::vector<std::string> tF1;
        tF1 = StrTokenize(TrimSpaces(tF[0]), '.');
        if((int)tF1.size() ==2)
        {
           ID tID = TrimSpaces(tF1[1]);
           if(std::find(allCifKeys["dataDesc"].begin(), allCifKeys["dataDesc"].end(),  
                     tID) != allCifKeys["dataDesc"].end())
           {
               dataDesc[tID] = tF[1];
           }
        }
        
    }
    
    void DictCifFile::getAtomInfo(std::vector<std::string> tF)
    {
        // check which entries exist in an atom
        if ((int)tF.size() ==1 && 
            tF[0].find("_chem_comp_atom.")!=std::string::npos)
        {
            std::vector<std::string> tF1;
            tF1 = StrTokenize(TrimSpaces(tF[0]), '.');
            
            if((int)tF1.size() ==2)
            {
                tF1[1] = TrimSpaces(tF1[1]);
                if(tF1[1].find("comp_id") !=std::string::npos)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("resName",curBlockLine));
                    //std::cout << curBlockLine << std::endl;
                    //curBlockLine++;
                }
                else if(tF1[1].find("atom_id") !=std::string::npos)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("id",curBlockLine));
                    //std::cout << curBlockLine << std::endl;
                    //curBlockLine++;
                }
                else if(tF1[1].find("type_symbol") !=std::string::npos)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("chemType",curBlockLine));
                    //std::cout << curBlockLine << std::endl;
                    //curBlockLine++;
                }
                else if(tF1[1].find("type_energy") !=std::string::npos)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("enerType",curBlockLine));
                    hasCCP4Type = true;
                    // std::cout << curBlockLine << std::endl;
                    //curBlockLine++;
                }
                else if(tF1[1].find("charge") !=std::string::npos 
                        && tF1[1].find("partial") == std::string::npos)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("charge",curBlockLine));
                    //std::cout << curBlockLine << std::endl;
                    //curBlockLine++;
                }
                else if (tF1[1].compare("x") ==0)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("x",curBlockLine) );
                    //curBlockLine++;
                    hasCoords = true;
                    //std::cout << curBlockLine << std::endl;
                }
                else if (tF1[1].compare("y") ==0)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("y",curBlockLine) );
                    // curBlockLine++;
                    //std::cout << curBlockLine << std::endl;
                }
                else if (tF1[1].compare("z") ==0)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("z",curBlockLine));
                    // curBlockLine++;
                    //std::cout << curBlockLine << std::endl;
                }
                else if (tF1[1].compare("model_Cartn_x") ==0)
                {
                    hasCoords = true;
                    hasProps["atom"].insert(std::pair<std::string, int>("model_Cartn_x",curBlockLine) );
                    //curBlockLine++;
                    //std::cout << curBlockLine << std::endl;
                }
                else if (tF1[1].compare("model_Cartn_y") ==0)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("model_Cartn_y", curBlockLine) );
                    // curBlockLine++;
                    //std::cout << curBlockLine << std::endl;
                }
                else if (tF1[1].compare("model_Cartn_z") ==0)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("model_Cartn_z",curBlockLine));
                    //curBlockLine++;
                    //std::cout << curBlockLine << std::endl;
                }
                else if (tF1[1].compare("pdbx_model_Cartn_x_ideal") ==0)
                {
                    hasCoords = true;
                    
                    hasProps["atom"].insert(std::pair<std::string, int>("x",curBlockLine) );
                    //curBlockLine++;
                    //std::cout << curBlockLine << std::endl;
                }
                else if (tF1[1].compare("pdbx_model_Cartn_y_ideal") ==0)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("y", curBlockLine) );
                    //curBlockLine++;
                    //std::cout << curBlockLine << std::endl;
                }
                else if (tF1[1].compare("pdbx_model_Cartn_z_ideal") ==0)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("z",curBlockLine));
                    
                    //curBlockLine++;
                    //std::cout << curBlockLine << std::endl;
                } 
                curBlockLine++;
            }
        }
        
        if ((int)tF.size() >2 && tF[0].find("#") ==std::string::npos)
        {
            itsCurAtom = new AtomDict();
            bool coordsDone = false;
            itsCurAtom->seriNum = itsCurAtomSeriNum;
            itsCurAtomSeriNum ++;
            //std::cout << "One new atom:" << std::endl;
            //std::cout << "Its seriNum  " <<  itsCurAtom->seriNum << std::endl;
            if (hasProps["atom"].find("resName") != hasProps["atom"].end())
            {
                itsCurAtom->existProps["resName"] =  hasProps["atom"]["resName"];
                itsCurAtom->resName = TrimSpaces(tF[hasProps["atom"]["resName"]]);
                //std::cout << "Its residue name: " << itsCurAtom->resName << std::endl;
            }
            if (hasProps["atom"].find("id") != hasProps["atom"].end())
            {
                itsCurAtom->existProps["id"] =  hasProps["atom"]["id"];
                cleanSymbol(tF[hasProps["atom"]["id"]], "\"");
                itsCurAtom->id = TrimSpaces(tF[hasProps["atom"]["id"]]);
                //std::cout << "Its ID: " << itsCurAtom->id << std::endl;
            }
            if (hasProps["atom"].find("chemType") != hasProps["atom"].end())
            {
                itsCurAtom->existProps["chemType"] =  hasProps["atom"]["chemType"];
                itsCurAtom->chemType = TrimSpaces(tF[hasProps["atom"]["chemType"]]);
                if (itsCurAtom->chemType.compare("D")==0)
                {
                    itsCurAtom->chemType="H";
                }
                
                //std::cout << "Its chemType : " << itsCurAtom->chemType << std::endl;
            }
            if (hasProps["atom"].find("enerType") != hasProps["atom"].end())
            {
                itsCurAtom->existProps["enerType"] = hasProps["atom"]["enerType"];
                itsCurAtom->ccp4Type = TrimSpaces(tF[itsCurAtom->existProps["enerType"]]);
                //std::cout << "Its energy type : " << itsCurAtom->enerType << std::endl;
            }
            if (hasProps["atom"].find("parCharge") != hasProps["atom"].end())
            {
                itsCurAtom->existProps["parChange"] = hasProps["atom"]["parCharge"];
                
                itsCurAtom->parCharge = StrToReal(tF[itsCurAtom->existProps["parChange"]]);
               // partial charge and formal charge are same in mmcif files in ccp4 monomer lib 
                // itsCurAtom->formalCharge = itsCurAtom->parCharge;
                //std::cout << "Its partialCharge :" 
                //        << itsCurAtom->parCharge
                //        << std::endl;
            }
            if (hasProps["atom"].find("charge") != hasProps["atom"].end())
            {
                itsCurAtom->existProps["charge"] = hasProps["atom"]["charge"];
                
                itsCurAtom->charge = StrToReal(tF[itsCurAtom->existProps["charge"]]);
                // partial charge and formal charge are same in mmcif files in ccp4 monomer lib 
                itsCurAtom->formalCharge  = itsCurAtom->charge;
                itsCurAtom->formalChargeI = itsCurAtom->charge;
                //std::cout << "Its formalCharge :" 
                //        << itsCurAtom->parCharge
                //        << std::endl;
            }
            if (hasProps["atom"].find("x") != hasProps["atom"].end())
            {
                
                itsCurAtom->existProps["x"] = hasProps["atom"]["x"];
                // std::cout << tF[itsCurAtom->existProps["x"]] << std::endl;
                std::string tSX = TrimSpaces(tF[itsCurAtom->existProps["x"]]);
                if (tSX.find("?") ==std::string::npos)
                {
                    itsCurAtom->coords[0]=StrToReal(tSX);
                    coordsDone = true;
                    // std::cout << "Its coord x : " << itsCurAtom->coords[0] << std::endl;
                }
            }
            if (hasProps["atom"].find("y") != hasProps["atom"].end())
            {
                
                itsCurAtom->existProps["y"] = hasProps["atom"]["y"];
                std::string tSY = TrimSpaces(tF[itsCurAtom->existProps["y"]]);
                if (tSY.find("?") ==std::string::npos)
                {
                    itsCurAtom->coords[1]=StrToReal(tSY);
                    // std::cout << "Its coord y : " << itsCurAtom->coords[1] << std::endl;
                }
            }
            if (hasProps["atom"].find("z") != hasProps["atom"].end())
            {
                itsCurAtom->existProps["z"] = hasProps["atom"]["z"];
                std::string tSZ = TrimSpaces(tF[itsCurAtom->existProps["z"]]);
                if (tSZ.find("?") ==std::string::npos)
                {
                    itsCurAtom->coords[2]=StrToReal(tSZ);
                    // std::cout << "Its coord z : " << itsCurAtom->coords[2] << std::endl;
                }    
            }
            if (!coordsDone)
            {
                if (hasProps["atom"].find("model_Cartn_x") != hasProps["atom"].end())
                {
                    itsCurAtom->existProps["x"] = hasProps["atom"]["model_Cartn_x"];
                    std::string tSX = TrimSpaces(tF[itsCurAtom->existProps["x"]]);
                    if (tSX.find("?") ==std::string::npos)
                    {
                        itsCurAtom->coords[0]=StrToReal(tSX);
                        // std::cout << "Its coord x : " << itsCurAtom->coords[0] << std::endl;
                        coordsDone = true;
                    }
                }
                if (hasProps["atom"].find("model_Cartn_y") != hasProps["atom"].end())
                {
                    itsCurAtom->existProps["y"] = hasProps["atom"]["model_Cartn_y"];
                    std::string tSY = TrimSpaces(tF[itsCurAtom->existProps["y"]]);
                    if (tSY.find("?") == std::string::npos)
                    {
                        itsCurAtom->coords[1]=StrToReal(tSY);
                        // std::cout << "Its coord y : " 
                        // << itsCurAtom->coords[1] << std::endl;
                    }
                }
                if (hasProps["atom"].find("model_Cartn_z") != hasProps["atom"].end())
                {
                    
                    itsCurAtom->existProps["z"] = hasProps["atom"]["model_Cartn_z"];
                    std::string tSZ = TrimSpaces(tF[itsCurAtom->existProps["z"]]);
                    if (tSZ.find("?") == std::string::npos)
                    {
                        
                        itsCurAtom->coords[2]=StrToReal(tSZ);
                        // std::cout << "Its coord z : " << itsCurAtom->coords[2] << std::endl;
                    }
                
                }
                // not have any kind of atom coords
                if (!coordsDone)
                {
                    // put a set of random coords to the atom if they do not exist
                    double r1 =  (double) rand()/RAND_MAX;
                    double r2 =  (double) rand()/RAND_MAX;
                    double r3 =  (double) rand()/RAND_MAX;
                    itsCurAtom->coords[0] = r1;
                    itsCurAtom->coords[1] = r2;
                    itsCurAtom->coords[2] = r3;
                }
            }
            
            allAtoms.push_back(*itsCurAtom);
            //std::cout << "Number of atoms now " << allAtoms.size() << std::endl;
            if (itsCurAtom->chemType == "H")
            {
                allHydroAtoms.push_back((int)allAtoms.size()-1);
            }
            
            delete itsCurAtom;
            itsCurAtom = NULL;
        }
       
    }
    
    void DictCifFile::setAtomsCCP4Type()
    {
        // should be done after rings and chiral centers are detected
        
        CCP4AtomType  aCPP4TypeTool(allAtoms, allRings);
        aCPP4TypeTool.setAllAtomsCCP4Type();
        for (int i=0; i < (int)aCPP4TypeTool.allAtoms.size(); i++)
        {
            allAtoms[i].ccp4Type = aCPP4TypeTool.allAtoms[i].ccp4Type;
            //std::cout << "Atom " << allAtoms[i].id 
            //          << " CCP4 atom type is " << allAtoms[i].ccp4Type 
            //          << std::endl;
        } 
    }
    
    int DictCifFile::atomPosition(ID tID)
    {
        for (int i=0; i<(int)allAtoms.size(); i++)
        {
            if(allAtoms[i].id.compare(tID)==0)
            {
                return i;
            }
        }
        
        return -1;
        
    }
    
    void DictCifFile::getBondAndConnect(std::vector<std::string> tF)
    {
        
        if ((int)tF.size() ==1 && 
            tF[0].find("_chem_comp_bond.")!=std::string::npos)
        {
            std::vector<std::string> tF1;
            tF1 = StrTokenize(TrimSpaces(tF[0]), '.');
            
            if((int)tF1.size() ==2)
            { 
                if(tF1[1].find("comp_id") !=std::string::npos)
                {
                    hasProps["bond"]["resName"] = curBlockLine;
                    curBlockLine++;
                }
                else if(tF1[1].find("atom_id_1") !=std::string::npos)
                {
                    hasProps["bond"]["atom1"] = curBlockLine;
                    curBlockLine++;
                }
                
                else if(tF1[1].find("atom_id_2") !=std::string::npos)
                {
                    hasProps["bond"]["atom2"] = curBlockLine;
                    curBlockLine++;
                }
                else if(tF1[1].find("type") !=std::string::npos 
                       || tF1[1].find("value_order") !=std::string::npos)
                {
                    hasProps["bond"]["order"] = curBlockLine;
                    curBlockLine++;
                }
                else if(TrimSpaces(tF1[1]).compare("value_dist") ==0)
                {
                    hasProps["bond"]["length"] = curBlockLine;
                    curBlockLine++;
                }
                else if (TrimSpaces(tF1[1]).compare("value_dist_esd")==0)
                {
                    hasProps["bond"]["sigLength"] = curBlockLine;
                    curBlockLine++;
                }
            }
        }
        
        //if ((int)tF.size() == (int)hasProps["bond"].size()
        //    && (int)tF.size() > 2 && tF[0].find("#") ==std::string::npos)
        if ((int)tF.size() > 2 && tF[0].find("#") ==std::string::npos)
        {
            itsCurBond = new BondDict();
            
            if (hasProps["bond"].find("resName") != hasProps["bond"].end())
            {
                itsCurBond->resName=tF[hasProps["bond"]["resName"]];
            }
            
            if (hasProps["bond"].find("atom1") != hasProps["bond"].end())
            {
                cleanSymbol(tF[hasProps["bond"]["atom1"]], "\"");
                std::string tId1 = TrimSpaces(tF[hasProps["bond"]["atom1"]]);
                itsCurBond->atoms.push_back(tId1);
            }
            
            if (hasProps["bond"].find("atom2") != hasProps["bond"].end())
            {
                cleanSymbol(tF[hasProps["bond"]["atom2"]], "\"");
                std::string tId2 = TrimSpaces(tF[hasProps["bond"]["atom2"]]);
                itsCurBond->atoms.push_back(tId2);   
            }
            
            if (hasProps["bond"].find("order") != hasProps["bond"].end())
            {
                itsCurBond->order   = tF[hasProps["bond"]["order"]];
                itsCurBond->orderNI = itsCurBond->order;
                itsCurBond->orderN = StrToOrder(itsCurBond->order);
            }
            
            if (hasProps["bond"].find("length") != hasProps["bond"].end())
            {
                itsCurBond->value = StrToReal(tF[hasProps["bond"]["length"]]);
            }
            
            if (hasProps["bond"].find("sigLength") != hasProps["bond"].end())
            {
                itsCurBond->sigValue = StrToReal(tF[hasProps["bond"]["sigLength"]]);
            }
            
            if ((int)itsCurBond->atoms.size() ==2)
            {
                itsCurBond->seriNum = itsCurBondSeriNum;
                itsCurBondSeriNum++; 
                // try using atom positions in vector allAtoms for the sake of calculation 
                int iPos1=0, iPos2=0;
                for (int i=0; i < (int) allAtoms.size(); i++ )   
                {
                    if(allAtoms[i].id.compare(itsCurBond->atoms[0])==0)
                    {
                        iPos1 = i;
                        itsCurBond->fullAtoms[allAtoms[i].id] = iPos1; 
                        itsCurBond->atomsIdx.push_back(iPos1);
                    }
                    else if(allAtoms[i].id.compare(itsCurBond->atoms[1])==0)
                    {
                        iPos2 = i;
                        itsCurBond->fullAtoms[allAtoms[i].id] = iPos2;
                        itsCurBond->atomsIdx.push_back(iPos2);
                    }
                }
            
                if (iPos1 >=0 && iPos1 < (int)allAtoms.size()
                        && iPos2 >=0 && iPos2 < (int)allAtoms.size())
                {
                    allAtoms[iPos1].connAtoms.push_back(iPos2);
                    if (std::find(allAtoms[iPos1].inBonds.begin(), 
                                  allAtoms[iPos1].inBonds.end(), itsCurBond->seriNum)
                            ==allAtoms[iPos1].inBonds.end())
                    {
                        allAtoms[iPos1].inBonds.push_back(itsCurBond->seriNum);
                    }
                    allAtoms[iPos2].connAtoms.push_back(iPos1);
                    if (std::find(allAtoms[iPos2].inBonds.begin(), 
                                  allAtoms[iPos2].inBonds.end(), itsCurBond->seriNum)
                            ==allAtoms[iPos2].inBonds.end())
                    {
                        allAtoms[iPos2].inBonds.push_back(itsCurBond->seriNum);
                    }
                }
            } 
            
            allBonds.push_back(*itsCurBond);
            delete itsCurBond;
            itsCurBond = NULL;
        }
    }
    
    void DictCifFile::setHydroAtomConnect()
    {
        for(std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA!=allAtoms.end(); iA++)
        {
            if(iA->chemType.compare("H")==0)
            {
                allHAtomIdx.push_back(iA->seriNum);
            }
            else
            {
                for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                    iNB !=iA->connAtoms.end(); iNB++)
                {
                    if (allAtoms[*iNB].chemType.compare("H")==0)
                    {
                        iA->connHAtoms.push_back(*iNB);
                    }
                }
            }
        }
        
        // Check
        /*
        for(std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA!=allAtoms.end(); iA++)
        {
            std::cout << "Atom " << iA->id << " connects to " 
                    << (int)iA->connHAtoms.size() << " H atoms " << std::endl;
            if ((int)iA->connHAtoms.size() !=0)
            {
                std::cout << "These atoms are : " << std::endl;
                for (std::vector<int>::iterator iH=iA->connHAtoms.begin();
                        iH !=iA->connHAtoms.end(); iH++)
                {
                    std::cout << "Atom " << allAtoms[*iH].id << std::endl;
                }
            }
        }
         */
       
        
    }
    void DictCifFile::setBondOrder()
    {
        // set the simplest cases
        setSingleOrder();
        
        if (hasConnect)
        {
            setBondOrderByBondList();
        }
        else if (hasCoords)
        {
            setBondOrderByCoords();
        }
        else
        {
            setBondOrderByChem();
        }
    }
    
    void DictCifFile::setSingleOrder()
    {
        // Set those bonds with certainty 
        // 1. atoms of 'H', 'F', 'CL', 'I', 'Br', 'B' have single covalent bonds
        std::map<std::string, int> covalenceMap;
        covalenceMap["H"]  = 1;
        covalenceMap["F"]  = 1;
        covalenceMap["CL"] = 1;
        covalenceMap["I"]  = 1;
        covalenceMap["BR"] = 1;
        covalenceMap["B"]  = 1;
        
        covalenceMap["D"]  = 1; // ?
        
        for (std::vector<BondDict>::iterator iB=allBonds.begin();
                iB !=allBonds.end(); iB++)
        {
            int idxA1 = iB->fullAtoms[iB->atoms[0]];
            int idxA2 = iB->fullAtoms[iB->atoms[1]];
            ID iD1 = allAtoms[idxA1].chemType;   
            ID iD2 = allAtoms[idxA2].chemType;
            // single bonds in the simplest cases
            std::map<std::string, int>::iterator iFind1;
            iFind1 = covalenceMap.find(iD1);
            if (iFind1 != covalenceMap.end() )
            {
                iB->order = 1;
            }
            else
            {
                std::map<std::string, int>::iterator iFind2;
                iFind2 = covalenceMap.find(iD2);
                if (iFind2 != covalenceMap.end())
                {
                    iB->order = 1;
                }
            }   
        }
    }
    
    void DictCifFile::setBondOrderByBondList()
    {
        std::map<ID,int> singleOrder;
        singleOrder["C"]  = 4;
        singleOrder["O"]  = 2;
        singleOrder["N"]  = 4;
        singleOrder["S"]  = 2;
        singleOrder["SE"] = 2;
        // allow other elements to be added in
        
        bool lAll = true;
        
        for (std::vector<BondDict>::iterator iB=allBonds.begin();
                iB !=allBonds.end(); iB++)
        {
            if (iB->order.empty())
            {
                int idxA1 = iB->fullAtoms[iB->atoms[0]];
                int idxA2 = iB->fullAtoms[iB->atoms[1]];
                ID iD1 = allAtoms[idxA1].chemType;   
                ID iD2 = allAtoms[idxA2].chemType;
                std::map<ID,int>::iterator iFind1 = singleOrder.find(iD1);        
                if (iFind1 !=singleOrder.end() && 
                    (int)allAtoms[idxA1].connAtoms.size() == singleOrder[iD1])
                {
                    iB->order = "single";
                }
                else 
                {
                    std::map<ID,int>::iterator iFind2 = singleOrder.find(iD2);
                    if (iFind2 !=singleOrder.end() &&
                        (int)allAtoms[idxA2].connAtoms.size() == singleOrder[iD2])
                    {
                        iB->order = "single";
                    }
                    else
                    {
                        lAll=false;
                    }
                }
            }
        }
        
        if (!lAll)
        {
            
        }
        
        
    }
    
    void DictCifFile::setBondOrderByChem()
    {
    }
    
    void DictCifFile::setBondOrderByCoords()
    {
        
    }
    
    void DictCifFile::addAtomSeriNumToBond()
    {
        for (int i=0; i < (int)allAtoms.size(); i++)
        {
            // add the seriNum to all bonds containing this angle.
            for (std::vector<BondDict>::iterator iBo=allBonds.begin();
                iBo !=allBonds.end(); iBo++)
            {
                if(allAtoms[i].id.compare(iBo->atoms[0])==0 ||
                        allAtoms[i].id.compare(iBo->atoms[1])==0)
                {
                    iBo->fullAtoms[allAtoms[i].id] = i;
                }
            }
        }
    }
    

    void DictCifFile::setAtomsCChemType()
    {
        
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA != allAtoms.end(); iA++)
        {
            // transfer the characters from 2nd position to low case
            if((int)iA->chemType.size()!=1)
            {
                ID tS = iA->chemType.substr(1);
                StrLower(tS);
                iA->chemType = iA->chemType[0] + tS;
            }
            
            // Embed planar info into the chemType to match COD definition 
            iA->cChemType = iA->chemType;
            if (iA->bondingIdx ==2)
            {
                 StrLower(iA->cChemType);
            }
            
            //std::cout << "Atom ID: " << iA->id << std::endl
            //        << "Atom chemType: " << iA->chemType << std::endl
            //        << "Atom chemType of COD classes: " << iA->cChemType << std::endl;
        }
    }
    
    void DictCifFile::setAtomsMetalType()
    {
        // setDefaultCoordGeos();
        
        ID metals[] = {"Li", "li", "Na", "na", "K",  "k",  "Rb", "rb", "Cs", "cs", "Fr", "fr",
                     "Be", "be", "Mg", "mg", "Ca", "ca", "Sr", "sr", "Ba", "ba", "Ra", "ra",
                     "Sc", "sc", "Y",  "y",
                     "Sb", "sb", "Te", "te", "Po", "po",
                     "Ti", "ti", "Zr", "zr", "Hf", "hf", "Rf", "rf",
                     "V",  "v"   "Nb", "nb", "Ta", "ta", "Db", "db", 
                     "Cr", "cr", "Mo", "mo", "W",  "w",  "Sg", "sg", 
                     "Mn", "mn", "Tc", "tc", "Re", "re", "Bh", "bh",  
                     "Fe", "fe", "Ru", "ru", "Os", "os", "Hs", "hs",   
                     "Co", "co", "Rh", "rh", "Ir", "ir", "Mt", "mt",  
                     "Ni", "ni", "Pd", "pd", "Pt", "pt", "Ds", "ds",  
                     "Cu", "cu", "Ag", "ag", "Au", "au", "Rg", "rg",   
                     "Zn", "zn", "Cd", "cd", "Hg", "hg",   
                     "Al", "al", "Ga", "ga", "In", "in", "Ti", "ti", 
                     "Sn", "sn", "Pb", "pb", "Bi", "bi"};
        // "Si", "si", "Ge", "ge", "As", "as",
        MetalTable.assign(metals, metals+115);
        /*
        std::cout << "Metal Elements :" << std::endl;
        for (std::vector<ID>::iterator iM =MetalTable.begin();
               iM !=MetalTable.end(); iM++)
        {
            std::cout << *iM << std::endl;
        }
        */
        
        
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA != allAtoms.end(); iA++)
        { 
              
            std::vector<ID>::iterator iFind = std::find(MetalTable.begin(), 
                                                        MetalTable.end(), iA->chemType);
            if (iFind !=MetalTable.end())
            {
                iA->isMetal = true;
                std::cout <<iA->id  << " of element " << iA->chemType 
                          << " is a metal atom." << std::endl;
                
                int cn = (int)iA->connAtoms.size();
                ID  ct = iA->chemType;
       
                if (iA->metalGeo.empty())
                {
                    std::map<ID, std::map<int,ID> >::iterator iFind1 
                                        =DefaultCoordGeos2.find(ct);
                    if (iFind1 !=DefaultCoordGeos2.end())
                    {
                        std::map<int,ID>::iterator iFind2
                        = DefaultCoordGeos2[ct].find(cn);
                        if (iFind2 !=DefaultCoordGeos2[ct].end())
                        {
                            iA->metalGeo=DefaultCoordGeos2[ct][cn];
                        }
                        else
                        {
                            iA->metalGeo = DefaultCoordGeos[cn];
                        }
                    }
                    else
                    {
                        iA->metalGeo = DefaultCoordGeos[cn];
                    }
                }
                
                std::cout << "Its coordination number is set to " << cn << std::endl
                          << "its default coordination geometry is set to " 
                          << iA->metalGeo << std::endl;
            }
        }
        
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA != allAtoms.end(); iA++)
        { 
            if (iA->isMetal)
            {
                std::cout << iA->id << " is a metal atom " << std::endl;
            }
            else
            {
                std::cout << iA->id << " is a organic atom " << std::endl;
            }
        }
        
    }
    
    void DictCifFile::addMissHydroAtoms()
    {
        // check valence and covalent bonds for each atom, add missing atom if 
        // necessary
        PeriodicTable   aPTab; 
        
        std::vector<ID> allMetals, pureMetals, metalloids;
        initMetalTab(pureMetals);
        initMetalloidTab(metalloids);
        for (std::vector<ID>::iterator iM=pureMetals.begin();
                iM !=pureMetals.end(); iM++)
        {
            allMetals.push_back(*iM);
        }
        for (std::vector<ID>::iterator iM=metalloids.begin();
                iM !=metalloids.end(); iM++)
        {
            allMetals.push_back(*iM);
        }
        
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA !=allAtoms.end(); iA++)
        {
            if (!isMetal(allMetals, iA->chemType))
            {
                // we have an organic atom, check its covalent bonds
                int nVBonds = 0, nMissH = 0;
            
                for (std::vector<int>::iterator iNB = iA->connAtoms.begin();
                        iNB !=iA->connAtoms.end(); iNB++)
                {
                    if(!isMetal(allMetals, allAtoms[*iNB].chemType))
                    {
                        nVBonds++;
                    }
                }
                // may need to consider case by case
                nMissH = aPTab.elements[iA->chemType]["val"] 
                        -aPTab.elements[iA->chemType]["pairedVE"]- nVBonds;
                
                if (nMissH)
                {
                    for (int iH =0; iH < nMissH; iH++)
                    {
                        // create a H atom and set its property 
                        AtomDict aHA;
                        aHA.id = "H" + iA->id;
                        aHA.seriNum = (int)allAtoms.size() + 1;
                        aHA.chemType = "H";
                        aHA.cChemType = "H";
                        aHA.connAtoms.push_back(atomPosition(iA->id));
                        aHA.resName = iA->resName;
                        aHA.isMetal = false;
                        
                        // add this atom to parent atom's connection, to allAtoms,
                        // and to allHydroAtoms
                        int iPos = (int)allAtoms.size();
                        iA->connAtoms.push_back(iPos);
                        allAtoms.push_back(aHA);
                        allHydroAtoms.push_back(iPos);
                    }
                } 
            }
            
        }
    }
    
    void DictCifFile::setAtomsVDWRadius()
    {
        PeriodicTable tPTable;
        
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA !=allAtoms.end(); iA++)
        {
            iA->radius= tPTable.elemProps[iA->chemType]["vdw"];
        }
    }
    
    void DictCifFile::setAtomsPartialCharges()
    {
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA !=allAtoms.end(); iA++)
        {
            
        }
    }
    
    void DictCifFile::CodAtomClassify(int dLev)
    {
        ringDetecting();
       
        
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA !=allAtoms.end(); iA++)
        {
            if ((int)iA->ringRep.size() >0)
            {
                std::cout << "\nAtom " << iA->id << " is in  " 
                        << (int)iA->ringRep.size() << " ring(s)" << std::endl;
                for (std::map<std::string, int>::iterator iRR=iA->ringRep.begin();
                        iRR != iA->ringRep.end(); iRR++)
                {
                    std::cout << "Ring: " << iRR->first << ", size : "
                            << iRR->second << std::endl;
                }
            }
        }
        
        // std::cout <<std::endl << "Output Atom COD classes now " << std::endl << std::endl;
        int iLev = 2;
        for (int i=0; i < (int)allAtoms.size(); i++)
        {
            setAtomCodClassName(allAtoms[i], allAtoms[i], iLev);
            std::cout <<std::endl << "For atom " << allAtoms[i].id << std::endl 
                      << "class is " << allAtoms[i].codClass << std::endl; 
                
        }
        /*
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
               iA !=allAtoms.end(); iA++)
        {
            setAtomCodClassName(*iA, *iA, iLev);
            std::cout <<std::endl << "For atom " << iA->id << std::endl 
                      << "class is " << iA->codClass << std::endl; 
                
        }
        */ 
    }
   
    void DictCifFile::detectPlaneGroups()
    {
        
        //groupOrgAtomsToPlanes();
        //groupMetAndLigandAtomsToPlanes();
   
        //Check
        std::cout<< "There are " << (int)allPlanes.size() 
                <<" Planes in the system" << std::endl;
        for (int i=0; i < (int)allPlanes.size(); i++)
        {
            std::cout<<"Plane " <<i+1 << " contains "
                    << (int)allPlanes[i].atoms.size() 
                    << " atoms. They are: "<< std::endl;
            for (std::map<ID, int>::iterator iAt=allPlanes[i].atoms.begin();
                    iAt!=allPlanes[i].atoms.end(); iAt++)
            {
                std::cout<<iAt->first << "\t";
            }
            std::cout<<std::endl;
        }
    }
    
    void DictCifFile::groupOrgAtomsToPlanes()
    {
        std::vector<PlaneDict> smalPls;
        
        // Find the smallest planes
        setSmallestPLs(smalPls);
        
        // Merge the above planes to the large plane groups
        mergeLargePLGroups(smalPls);
        
    }
    
    void DictCifFile::setSmallestPLs(std::vector<PlaneDict>& tSmaPls)
    {
        
        for (int i=0; i < (int)allAtoms.size(); i++)
        {
            if (allAtoms[i].bondingIdx == 2)
            {
                PlaneDict tSmaPl;
                ID        aID     = allAtoms[i].id;
                tSmaPl.archID  = aID;
                tSmaPl.archPos = i;
                tSmaPl.atoms[aID] = i;
                
                for (int j=0; j< (int)allAtoms[i].connAtoms.size(); j++)
                {
                    int tPos  = allAtoms[i].connAtoms[j];
                    if (!allAtoms[tPos].isMetal)
                    {
                        ID  tId   = allAtoms[tPos].id;
                        tSmaPl.atoms[tId] = tPos;
                    }
                }
                tSmaPls.push_back(tSmaPl);
            }
        }
        
        std::cout <<"There are " << (int)tSmaPls.size() 
                  << " sp2 planes " << std::endl;
        for (std::vector<PlaneDict>::iterator iP=tSmaPls.begin();
                iP !=tSmaPls.end(); iP++)
        {
            std::cout <<"small plane: " << std::endl;
            std::cout << "Center atom " << iP->archID << std::endl;
            for (std::map<ID, int>::iterator iA=iP->atoms.begin();
                    iA !=iP->atoms.end(); iA++)
            {
                if (iA->first != iP->archID )
                {
                    std::cout << "Other atom " << iA->first << std::endl;
                }
            }
        }
        
        
    }
    
    void DictCifFile::mergeLargePLGroups(std::vector<PlaneDict>& tSmaPls)
    {
        //std::cout <<"There are " << (int)tSmaPls.size() 
        //          << " sp2 planes " << std::endl;
        std::map<int, std::vector<int> > allPlsClasses;
        
        // assign a e-id to each plane set
        // allPlsClasses[e-classId] = vector of small plane indexes 
        for (int i=0; i < (int)tSmaPls.size(); i++)
        {
            for (std::map<ID,int>::iterator iSA=tSmaPls[i].atoms.begin();
                   iSA!=tSmaPls[i].atoms.end(); iSA++ )
            {
                allPlsClasses[i].push_back(iSA->second);
            }
        }
             
        // get equiv-class by ring-relation
        
        for (int i=0; i < (int)tSmaPls.size(); i++)
        {
            bool iMerge = false;
            /*
            std::cout << "small plane " << i << std::endl;
            std::cout << "Archor atoms: " <<tSmaPls[i].archID << std::endl;
            for (std::map<ID, int>::iterator iSA1=tSmaPls[i].atoms.begin();
                   iSA1 != tSmaPls[i].atoms.end(); iSA1++)
            {
                std::cout << iSA1->first << "\t";
            }
            std::cout<<std::endl;
            */
            
            for (int j=i+1; j <(int)tSmaPls.size(); j++)
            {
                /*
                std::cout << "linked small plane " << j << std::endl;
                std::cout << "Archor atoms: " << tSmaPls[j].archID << std::endl;
                for (std::map<ID, int>::iterator iSA2=tSmaPls[j].atoms.begin();
                   iSA2 != tSmaPls[j].atoms.end(); iSA2++)
                {
                    std::cout << iSA2->first << "\t";
                }
                std::cout<<std::endl;
                
                 */
                
                std::map<int, std::vector<int> >::iterator tFindM;
                tFindM = allPlsClasses.find(i);
                
                if (tFindM != allPlsClasses.end())
                {
                    //std::cout << "test ring-related " << std::endl;
                    //std::cout << isInSameRing(tSmaPls[i],tSmaPls[j]) << std::endl;
                   
                    if(isInSameRing(tSmaPls[i],tSmaPls[j]))
                    {   
                        // std::cout << "Merged planes " << i <<"  "<< j << std::endl;
                        for (std::vector<int>::iterator iI = allPlsClasses[i].begin();
                                iI !=allPlsClasses[i].end(); iI++ )
                        {
                            std::vector<int>::iterator tFindV;
                            tFindV = std::find(allPlsClasses[j].begin(), 
                            allPlsClasses[j].end(), *iI);
                            if (tFindV ==allPlsClasses[j].end())
                            {
                                allPlsClasses[j].push_back(*iI); 
                            }
                        }
                        //std::cout << "Plane " << j << "now has atoms " 
                        //          << (int)allPlsClasses[j].size() << std::endl;
                        iMerge = true;        
                    }    
                }
            }
            
            if (iMerge)
            {
                allPlsClasses.erase(i);
            }
        }  
        
        // Further
        
        std::vector<int> delKeys;
        for(std::map<int, std::vector<int> >::iterator iPl1=allPlsClasses.begin();
                iPl1!=allPlsClasses.end(); iPl1++)
        {
            for(std::map<int, std::vector<int> >::iterator iPl2=allPlsClasses.begin();
                iPl2!=allPlsClasses.end(); iPl2++)
            {
                if (iPl1->first < iPl2->first)
                {   
                    if(furtherM(iPl1->second, iPl2->second))
                    {    
                        
                        if((int)iPl1->second.size() >= (int)iPl2->second.size())
                        {
                            for (std::vector<int>::iterator iV=iPl2->second.begin();
                                    iV!=iPl2->second.end(); iV++)
                            {
                                
                                std::vector<int>::iterator tFindV;
                                tFindV = std::find(iPl1->second.begin(), 
                                         iPl1->second.end(), *iV);
                                if (tFindV ==iPl1->second.end())
                                {
                                    iPl1->second.push_back(*iV); 
                                }
                            }
                            delKeys.push_back(iPl2->first);   
                        }
                        else
                        {
                            for (std::vector<int>::iterator iV=iPl1->second.begin();
                                    iV!=iPl1->second.end(); iV++)
                            {
                                std::vector<int>::iterator tFindV;
                                tFindV = std::find(iPl2->second.begin(), 
                                         iPl2->second.end(), *iV);
                                if (tFindV ==iPl2->second.end())
                                {
                                    iPl2->second.push_back(*iV); 
                                }
                            }
                            delKeys.push_back(iPl1->first);
                        }
                    }
                }
            }
        }
        
        for(std::vector<int>::iterator iDel=delKeys.begin();
                iDel!=delKeys.end(); iDel++)
        {
            if(allPlsClasses.find(*iDel) != allPlsClasses.end())
            {
                allPlsClasses.erase(*iDel);
            }
        }
       
        
        // Put equiv-classes in a vector of planeDicts 
        for (std::map<int, std::vector<int> >::iterator iPC=allPlsClasses.begin();
                iPC!=allPlsClasses.end(); iPC++ )
        {
            PlaneDict aPl;
            
            for (std::vector<int>::iterator iId=iPC->second.begin();
                   iId !=iPC->second.end(); iId++ )
            {
                aPl.atoms[allAtoms[*iId].id]=*iId;
            }  
            allPlanes.push_back(aPl);
        } 
    }

    bool DictCifFile::isInSameRing(PlaneDict & tP1, PlaneDict & tP2)
    {
        if (std::find(allAtoms[tP1.archPos].connAtoms.begin(),
                allAtoms[tP1.archPos].connAtoms.end(), tP2.archPos)
                ==allAtoms[tP1.archPos].connAtoms.end())
        {
            // std::cout << "Not in " << tP1.archID << "Neighbor " << std::endl; 
            return false;
        } 
        
        if ((int)allAtoms[tP1.archPos].ringRep.size())
        {
                for (std::map<std::string, int>::iterator iM1=allAtoms[tP1.archPos].ringRep.begin();
                        iM1 !=allAtoms[tP1.archPos].ringRep.end(); iM1++)
                {
                    // std::cout<< "1: " << iM1->first << std::endl;
                    if ((int)allAtoms[tP2.archPos].ringRep.size())
                    {
                       /* 
                        for (std::map<std::string, int>::iterator iM2=allAtoms[tP2.archPos].ringRep.begin();
                                 iM2 !=allAtoms[tP2.archPos].ringRep.end(); iM2++)
                        {
                            std::cout<< "2: " << iM2->first << std::endl;
                        }
                        */
                            
                        if(allAtoms[tP2.archPos].ringRep.find(iM1->first) !=
                                   allAtoms[tP2.archPos].ringRep.end())
                        {
                            if (allRings[iM1->first][0].isPlanar)
                            {
                                return true;
                            }
                        }
                    }
                }
            
        }
        
        return false;   
    }
    
    bool DictCifFile::furtherM(std::vector<int> &tV1, std::vector<int> &tV2)
    {
        int nFind = 0;
        for(std::vector<int>::iterator iV=tV1.begin();
                iV !=tV1.end(); iV++)
        {
           std::vector<int>::iterator tFindV;
           tFindV = std::find(tV2.begin(), tV2.end(), *iV);
           if (tFindV !=tV2.end())
           {
               nFind++; 
           } 
        }
        
        if (nFind >=3)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
        
    // Ring related 
    void DictCifFile::ringDetecting()
    {
        // tempo, set ring size to 6
        int maxSize = 7;
        // std::vector<AtomDict> atomsInPath;
        std::map<int, ID>  atomsInPath;
        std::map<int, ID>  atomsSeen;
       
        // 1. loops beginning from all atoms in the system
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA !=allAtoms.end(); iA++)
        {
   
            //std::cout << "-----------------" << std::endl;
            //std::cout << "starting atom  " << iA->id << std::endl;
            //std::cout << "-----------------" << std::endl;
            
            //if((int)tempIDs.size() !=0)
            //{
            //    tempIDs.clear();
            //}
            
            int preSeriNum = -999;
            int startLev   = 1;
            atomsInPath.clear();
            atomsSeen.clear();
            // atomsInPath.push_back((*iA));
            // atomsInPath.insert(std::pair<int, ID>(iA->seriNum, iA->chemType))
            // tempIDs.insert(std::pair<int, ID>(iA->seriNum, iA->id));
            
            // 2. loop from its bonded atoms recursively
            // checkOnePathSec(atomsInPath, tempIDs, *iA, maxSize, iA);
            // checkOnePathSec(atomsInPath, *iA, maxSize, iA);
            checkOnePathSec(*iA, maxSize, iA, preSeriNum,  startLev, atomsSeen, atomsInPath);
            
            /*
            if ((int)iA->ringRep.size() >0)
            {
                
                std::cout << "\nAtom " << iA->id << " is in "
                        << (int)iA->ringRep.size() << " rings "
                        << std::endl;
                for (std::map<std::string, int>::iterator iMa = iA->ringRep.begin();
                       iMa != iA->ringRep.end(); iMa++ )
                {
                    std::cout << "Ring: " << iMa->first << "; Size " << iMa->second 
                            <<std::endl;
                } 
            } 
            else
            {
                std::cout << "\nAtom " << iA->id << " is in no ring" <<std::endl;
            }
            */
           
            // std::cout << "finish atom  " << iA->seriNum << std::endl;
            /*
            if ((int)atomsInPath.size() >0)
            {
                for (std::vector<Ring>::iterator iR=atomsInPath[0].inRings.begin();
                      iR != atomsInPath[0].inRings.end(); iR++)
                {
                    iA->inRings.push_back(*iR);
                }
            }
             */
            
        
        }
        
        /*
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA !=allAtoms.end(); iA++)
        {
            if ((int)iA->ringRep.size() >0)
            {
                std::cout << "\nAtom " << iA->id << " is in  " 
                        << (int)iA->ringRep.size() << " ring(s)" << std::endl;
                for (std::map<std::string, int>::iterator iRR=iA->ringRep.begin();
                        iRR != iA->ringRep.end(); iRR++)
                {
                    std::cout << "Ring: " << iRR->first << ", size : "
                            << iRR->second << std::endl;
                }
            }
        }
         */
    }
    
   
    // Version 2 
    void DictCifFile::checkOnePathSec(AtomDict                & curAto,
                                      int                       iMax,
                                      std::vector<AtomDict>::iterator iOriAto,
                                      int                       SeriNumPreAto,  
                                      int                       curLev,
                                      std::map<int, ID>       & seenAtomIDs,
                                      std::map<int, ID>       & atomIDsInPath)
    {
       
        if ( curLev <iMax )
        {  
            int NachbarpunkteDetected = 0;
            
            // Check Nachbarpunkte
            for (std::vector<int>::iterator tNBA=curAto.connAtoms.begin();
                    tNBA != curAto.connAtoms.end(); tNBA++)
            {
                int tSeriNum = allAtoms[*tNBA].seriNum;
                
                // Find  Nachbarpunkte of the current path if the current atom
                // (1) should not be the the atom beginning the path. (it is a ring)
                // (2) should not be the one just walked through in the last step
                // (3) is in a list of atoms we have seen
                if (tSeriNum != iOriAto->seriNum && tSeriNum != SeriNumPreAto 
                        && seenAtomIDs.find(tSeriNum) !=seenAtomIDs.end())
                {
                    
                    NachbarpunkteDetected = 1;
                    
                    // found Nachbarpunkte, stop this path
                    /*
                    
                    std::cout << "atom : " <<  allAtoms[*tNBA].id << std::endl;
                    
                    std::cout << "Nachbarpunkte found, stop this path " << std::endl;
                    for (std::map<int, ID>::iterator iS=seenAtomIDs.begin();
                            iS != seenAtomIDs.end(); iS++)
                    {
                        std::cout << "atom : " << iS->first 
                                << " : " << iS ->second << std::endl;
                    }
                    */              
                }  
            }
            
            // check if a ring is formed 
            if (!NachbarpunkteDetected)
            {
                for (std::vector<int>::iterator tNBA=curAto.connAtoms.begin();
                    tNBA != curAto.connAtoms.end(); tNBA++)
                {
                    int tSeriNum = allAtoms[*tNBA].seriNum;
                    
                    if (tSeriNum == iOriAto->seriNum && tSeriNum != SeriNumPreAto
                        && curLev > 2)
                    {
                        //std::cout << iOriAto->id << " : "
                        //          << curAto.id << " : " << allAtoms[*tNBA].id
                        //          << " : " << SeriNumPreAto 
                        //          << " Find a ring." << std::endl;    
                        //   sort the atoms in the seenAtom vector and 
                        //   check if this ring is already found (the same 
                        //   ring of the same atom should be found twice at
                        //   least because of the walking algorithm.
                        // FIND A RING !
                        atomIDsInPath.insert(std::pair<int,ID>(curAto.seriNum,curAto.id));
                        std::list<std::string> tAllIds;
                        std::vector <AtomDict> ttAtoms;
                        for (std::map<int, ID>::iterator iSee = atomIDsInPath.begin();
                                iSee != atomIDsInPath.end(); iSee++)
                        {
                            //if (iSee->first != iOriAto->seriNum)
                            //{
                            tAllIds.push_back(iSee->second);
                            int posIdx = atomPosition(iSee->second);
                            ttAtoms.push_back(allAtoms[posIdx]);
                                // tRing.atoms.push_back(*tAt);
                                //std::cout << iSee->second << std::endl;   
                            //}
                        }
                        RingDict aRingDict(ttAtoms);                 
                        
                        tAllIds.sort(compareNoCase);
                        //std::string tRepStr(iOriAto->id);
                        std::string tRepStr;
                        for (std::list<std::string>::iterator iAI =tAllIds.begin();
                                    iAI != tAllIds.end(); iAI++)
                        {
                            tRepStr.append(*iAI);
                        }
                            
                        iOriAto->ringRep[tRepStr] = (int)atomIDsInPath.size();
                        
                        aRingDict.rep = tRepStr;
                        std::map<ID, std::vector<RingDict> >::iterator iFindRing=allRings.find(tRepStr);
                        if (iFindRing == allRings.end())
                        {
                            allRings[tRepStr].push_back(aRingDict);
                        }
                        
                        atomIDsInPath.erase(curAto.seriNum);  
                        NachbarpunkteDetected = 1;
                        
                    }
                    
                }
            }
              
            if (! NachbarpunkteDetected)
            {
                // found no Nachbarpunkte and no ring
                // descend the new path in the neighborhood graph:
                
                
                int tNewLev = curLev + 1;
                seenAtomIDs.insert(std::pair<int,ID>(curAto.seriNum,curAto.id));
                if (tNewLev < iMax)
                {
                   /*
                        std::cout << "atom " << curAto.id 
                                  << " finds no Nachbarpunkte in it neighbor  " 
                                  << std::endl << "Descend into the next atom "
                                  << std::endl;
                    
                    */
                    
                    atomIDsInPath.insert(std::pair<int,ID>(curAto.seriNum,curAto.id));
                    for (std::vector<int>::iterator tNBA=curAto.connAtoms.begin();
                         tNBA != curAto.connAtoms.end(); tNBA++)
                    {
                        if(curLev==1)
                        {
                            // tempo list of atoms in a path
                            if((int)seenAtomIDs.size() !=0)
                            {
                                seenAtomIDs.clear();
                            }
                            if((int)atomIDsInPath.size() !=0)
                            {
                                atomIDsInPath.clear();
                            }
                            //std::cout << "after clear, the size is " 
                            //          << (int)seenAtomIDs.size() << std::endl;
                            seenAtomIDs.insert(std::pair<int,ID>(curAto.seriNum,curAto.id));
                            atomIDsInPath.insert(std::pair<int,ID>(curAto.seriNum,curAto.id));
                        }
                        if (SeriNumPreAto != allAtoms[*tNBA].seriNum)
                        {
                            /*
                                std::cout << std::endl << "Current size " << curLev << std::endl;
                                std::cout << "Orig atom : " << iOriAto->id
                                          << " Prev atom : " << SeriNumPreAto
                                          << " Curent atom :  " << curAto.id 
                                          << std::endl << std::endl;
                                std::cout << "NB atom : " << allAtoms[*tNBA].id << std::endl;  
                            */
                            // This is a new atom and append the atom to seen-atom list
                            // and call function checkOnePathSec() recursively
                            int tPreSeriNum = curAto.seriNum; 
                            checkOnePathSec(allAtoms[*tNBA], iMax, iOriAto, tPreSeriNum, 
                                        tNewLev, seenAtomIDs, atomIDsInPath);
                        }
                    }
                    atomIDsInPath.erase(curAto.seriNum);
                    seenAtomIDs.erase(curAto.seriNum);
                }
                atomIDsInPath.erase(curAto.seriNum);
                seenAtomIDs.erase(curAto.seriNum);
            }
        }
        else
        {
            atomIDsInPath.erase(curAto.seriNum);
            seenAtomIDs.erase(curAto.seriNum);
        }
    }
    
    
    void DictCifFile::setAtomCodClassName(AtomDict & tAtom,
                                          AtomDict & tOriAtom,
                                          int tLev)
    {
        
        if (tLev==1)
        {
            tAtom.codClass = "";
            tAtom.codClass.append(tAtom.chemType);
            outRingSec(tAtom);
            
            std::string tStr;
            std::list<std::string> tStrList;
            std::map<ID, int> comps;
            
            //tStrList.push_back(tAtom.chemType);
            
            //std::cout << "Id list size " << (int) tStrList.size() << std::endl;
            // just get immediate neighbor atom ID
            for (std::vector<int>::iterator tNBAtom=tAtom.connAtoms.begin();
                    tNBAtom != tAtom.connAtoms.end(); tNBAtom++)
            {
                if(allAtoms[*tNBAtom].seriNum != tOriAtom.seriNum)
                {
                    // tStrList.push_back(allAtoms[*tNBAtom].chemType);
                    if(comps.find(allAtoms[*tNBAtom].chemType) != comps.end())
                    {
                        comps[allAtoms[*tNBAtom].chemType] += 1;
                    }
                    else
                    {
                        comps[allAtoms[*tNBAtom].chemType] = 1; 
                    }
                }
            }
            
            for (std::map<ID, int>::iterator iMa=comps.begin();
                    iMa !=comps.end(); iMa++)
            {
                std::string s1, s2;
                s1 = iMa->first + IntToStr(iMa->second);
                for (int i=0; i < iMa->second; i++)
                {
                    s2.append(iMa->first);
                }
                if ((int)s1.size() < (int)s2.size())
                {
                    tStrList.push_back(s1);
                }
                else
                {
                    tStrList.push_back(s2);
                }
            }
            
            tStrList.sort(compareNoCase);
            // std::cout << "sort Id list size " << (int) tStrList.size() << std::endl;
            for (std::list<std::string>::iterator iL = tStrList.begin();
                    iL != tStrList.end(); iL++)
            {
                tStr.append(*iL);
            }
            
            //std::cout << "the final str size " << (int) tStr.size() << std::endl;
            tAtom.nbRep.push_back(tStr);
            tAtom.codClass.append(tStr);
        }
        else if(tLev==2)
        {
            tAtom.codClass = "";
            tAtom.codClass.append(tAtom.chemType);
            outRingSec(tAtom);
            //std::cout << "Atom " << tAtom.id << " its COD ring section " 
            //        <<  tAtom.codClass << std::endl;
            
            int lowLev = tLev - 1;
            std::map<std::string, int> tIdMap;
            for (std::vector<int>::iterator tNBA=tAtom.connAtoms.begin();
                    tNBA != tAtom.connAtoms.end(); tNBA++)
            {
                setAtomCodClassName(allAtoms[*tNBA], tOriAtom, lowLev);
                /*
                std::list<std::string> tStrList;
                std::string tStr(allAtoms[*tNBA].chemType);
                tStr.append(outRingSecStr(allAtoms[*tNBA]));
                
                for (std::vector<int>::iterator tNNBA=allAtoms[*tNBA].connAtoms.begin();
                        tNNBA != allAtoms[*tNBA].connAtoms.end(); tNNBA++)
                {
                    if(allAtoms[*tNNBA].id.compare(tAtom.id) !=0)
                    {
                        tStrList.push_front(allAtoms[*tNNBA].chemType);
                    }
                }
                tStrList.sort(compareNoCase);
                for (std::list<std::string>::iterator iSL=tStrList.begin();
                        iSL != tStrList.end(); iSL++)
                {
                    tStr.append(*iSL);
                }
                */
                
                if(tIdMap.find(allAtoms[*tNBA].codClass) !=tIdMap.end())
                {
                    tIdMap[allAtoms[*tNBA].codClass]++;
                    
                }
                else
                {
                    tIdMap[allAtoms[*tNBA].codClass] = 1;
                }
                 
            }
            
            sortMap  tSMap;
            std::vector<sortMap> tVec;
            
            for (std::map<std::string, int>::iterator tM=tIdMap.begin();
                   tM !=tIdMap.end(); tM++)
            {
                tSMap.key = tM->first;
                tSMap.val = tM->second;
                tVec.push_back(tSMap);
            }
            
            std::sort(tVec.begin(),tVec.end(), desSortMapKey);
            
            // check
            /*
            if (tAtom.id == "B4")
            {
               std::cout << "After sorting " << std::endl;
               for (std::vector<sortMap>::iterator iV=tVec.begin();
                    iV != tVec.end(); iV++)
               {
                    std::cout << " key: " << iV->key << " value : "
                          << iV->val << std::endl;
               }
            }
            */
            for(std::vector<sortMap>::iterator iV=tVec.begin();
                    iV !=tVec.end(); iV++)
            {
                
                if(iV->val ==1 && (int)iV->key.length()==1)
                {
                    tAtom.codClass.append(iV->key);
                }
                else if (iV->val ==1)
                {
                    tAtom.codClass.append("("+iV->key+")");
                }
                else
                {
                    tAtom.codClass.append("(" + iV->key + ")" + IntToStr(iV->val));
                }
            }
            
            //std::cout<<"For atom " << tAtom.id << " : " << std::endl;
            //std::cout << "Its COD class is : " << tAtom.codClass 
            //          << std::endl <<std::endl;
        }
    }
    
    
    void DictCifFile::outRingSec(AtomDict &tAtom)
    {
        int numRings = (int)tAtom.ringRep.size();
       
        if (numRings)
        {
            
            std::map<int, int> sizeMap;
            
            
            for (std::map<std::string, int>::iterator iMR=tAtom.ringRep.begin();
                    iMR != tAtom.ringRep.end(); iMR++)
            {
                if (sizeMap.find(iMR->second) ==sizeMap.end())
                {
                    sizeMap[iMR->second] =1;
                }
                else
                {
                    sizeMap[iMR->second]++;
                }   
            }
            
            tAtom.codClass.append("[");
            int i =0;
            int j = (int)sizeMap.size();
            for (std::map<int, int>::iterator iSMa=sizeMap.begin();
                    iSMa != sizeMap.end(); iSMa++)
            {
                std::string tSize = IntToStr(iSMa->first);
                std::string tNum  = IntToStr(iSMa->second);
                
                if(iSMa->second >= 3)
                {
                    tAtom.codClass.append(tNum + "x" + tSize);
                }
                else if (iSMa->second==2)
                {
                    tAtom.codClass.append( tSize + "," + tSize);
                }   
                else if (iSMa->second==1)
                {
                    tAtom.codClass.append(tSize);
                }
                       
                
                if(i != j-1)
                {
                    tAtom.codClass.append(",");
                }
                else
                {
                    tAtom.codClass.append("]");
                }
                    
                i++;
            }
        }
    }
    
    std::string  DictCifFile::outRingSecStr(AtomDict &tAtom)
    {
        std::string tS1 = "";
        int numRings = (int)tAtom.ringRep.size();
        
        if (numRings)
        {
            
            std::map<int, int> sizeMap;
            
            
            for (std::map<std::string, int>::iterator iMR=tAtom.ringRep.begin();
                    iMR != tAtom.ringRep.end(); iMR++)
            {
                if (sizeMap.find(iMR->second) ==sizeMap.end())
                {
                    sizeMap[iMR->second] =1;
                }
                else
                {
                    sizeMap[iMR->second]++;
                }   
            }
            
            tS1.append("[");
            int i =0;
            int j = (int)sizeMap.size();
            for (std::map<int, int>::iterator iSMa=sizeMap.begin();
                    iSMa != sizeMap.end(); iSMa++)
            {
                std::string tSize = IntToStr(iSMa->first);
                std::string tNum  = IntToStr(iSMa->second);
               
                if(iSMa->second !=1)
                {
                    tS1.append(tNum + "x" + tSize);
                }
                else
                {
                    tS1.append(tSize);
                }
                
                if(i != j-1)
                {
                    tS1.append(",");
                }
                else
                {
                    tS1.append("]");
                }
                    
                i++;
            }
        }
        
        return tS1;
    }
      
    void DictCifFile::getAngleInfo(std::vector<std::string> tF)
    {   
        
        if ((int)tF.size() ==1 && 
            tF[0].find("_chem_comp_angle.")!=std::string::npos)
        {
            std::vector<std::string> tF1;
            tF1 = StrTokenize(TrimSpaces(tF[0]), '.');
            if((int)tF1.size() ==2)
            {
                if(tF1[1].find("comp_id") !=std::string::npos)
                {
                    hasProps["angle"].insert(std::pair<std::string, int>("resName",curBlockLine));
                    // std::cout << curBlockLine << std::endl;
                    curBlockLine++;
                }
                else if(tF1[1].find("atom_id_1") !=std::string::npos)
                {
                    hasProps["angle"].insert(std::pair<std::string, int>("atom_id_1",curBlockLine));
                    //std::cout << curBlockLine << std::endl;
                    curBlockLine++;
                }
                else if(tF1[1].find("atom_id_2") !=std::string::npos)
                {
                    hasProps["angle"].insert(std::pair<std::string, int>("atom_id_2",curBlockLine));
                    //std::cout << curBlockLine << std::endl;
                    curBlockLine++;
                }
                else if(tF1[1].find("atom_id_3") !=std::string::npos)
                {
                    hasProps["angle"].insert(std::pair<std::string, int>("atom_id_3",curBlockLine));
                    // std::cout << curBlockLine << std::endl;
                    curBlockLine++;
                }
                else if(tF1[1].find("value_angle") !=std::string::npos)
                {
                    hasProps["angle"].insert(std::pair<std::string, int>("value", curBlockLine));
                    //std::cout << curBlockLine << std::endl;
                    curBlockLine++;
                }
                else if(tF1[1].find("value_angle_esd") !=std::string::npos)
                {
                    hasProps["angle"].insert(std::pair<std::string, int>("value_esd",curBlockLine));
                    //std::cout << curBlockLine << std::endl;
                    curBlockLine++;
                }  
                
            }
            
        }
        
        //  remember the residue name is ignored so we have
        if ((int)tF.size() == (int)hasProps["angle"].size()+1 
                && (int)tF.size() >2 && tF[0].find("#") ==std::string::npos)
        {
            itsCurAngle = new AngleDict();
            itsCurAngle->seriNum = itsCurAngleSeriNum;
            itsCurAngleSeriNum++;
            
            
            
            //  No need for that part, could rewrite them
            //if (hasProps["angle"].find("resName") != hasProps["angle"].end())
            //{
            //    itsCurAngle->resName=tF[hasProps["angle"]["resName"]];
            //}
            
            
            if (hasProps["angle"].find("atom_id_1") != hasProps["angle"].end())
            {
                cleanSymbol(tF[hasProps["angle"]["atom_id_1"]], "\"");
                itsCurAngle->atoms.push_back(StrToInt(TrimSpaces(tF[hasProps["angle"]["atom_id_1"]])));
            }
            
            if (hasProps["angle"].find("atom_id_2") != hasProps["angle"].end())
            {
                cleanSymbol(tF[hasProps["angle"]["atom_id_2"]], "\"");
                itsCurAngle->atoms.push_back(StrToInt(TrimSpaces(tF[hasProps["angle"]["atom_id_2"]])));
            }
            
            if (hasProps["angle"].find("atom_id_3") != hasProps["angle"].end())
            {
                cleanSymbol(tF[hasProps["angle"]["atom_id_3"]], "\"");
                itsCurAngle->atoms.push_back(StrToInt(TrimSpaces(tF[hasProps["angle"]["atom_id_3"]])));
            }
            
            // temp keep CCP4 dictionary values first 
            if (hasProps["angle"].find("value") != hasProps["angle"].end())
            {
                itsCurAngle->value = StrToReal(tF[hasProps["angle"]["value"]]);
            }
            
            if (hasProps["angle"].find("value_esd") != hasProps["angle"].end())
            {
                itsCurAngle->sigValue = StrToReal(tF[hasProps["bond"]["value_esd"]]);
            }
            
            allAngles.push_back(*itsCurAngle);
            delete itsCurAngle;
            itsCurAngle = NULL;
        } 
        
    }
    
   // setAllAngles() may not needed in the future
    void DictCifFile::setAllAngles()
    {
        
        for(int i=0; i < (int)allAtoms.size(); i++)
        {   
            for (int j=0; j < (int)allAtoms[i].connAtoms.size(); j++)
            {
                for (int k=j+1; k < (int)allAtoms[i].connAtoms.size(); k++)
                {
                    int i1 = allAtoms[i].connAtoms[j];
                    int i2 = allAtoms[i].connAtoms[k];
                    //std::cout << "Angle between " << allAtoms[i].id 
                    //          << "(center) and " << allAtoms[i1].id
                    //          << " and " << allAtoms[i2].id << std::endl;
                        
                    AngleDict aAng;
                    std::vector<int> tVec;
                    aAng.anchorID  = allAtoms[i].id;
                    aAng.anchorPos = i;
                    aAng.atoms.push_back(i);
                       
                    if ((int) allAtoms[i1].connAtoms.size() >=
                        (int) allAtoms[i1].connAtoms.size())
                    {
                        aAng.atoms.push_back(i1);
                        aAng.atoms.push_back(i2);
                        tVec.push_back(i1);
                        tVec.push_back(i2);
                    }
                    else
                    {
                        aAng.atoms.push_back(i2);
                        aAng.atoms.push_back(i1);
                        tVec.push_back(i1);
                        tVec.push_back(i2);
                    }
                        
                    aAng.value        = 0.0;
                    aAng.sigValue     = 3.0;
                    aAng.numCodValues = 0;
                    allAngles.push_back(aAng);
                    allAnglesIdxs[i].push_back(tVec);
                }
            }
        }
        
        // Check 
        std::cout << "There are " << (int)allAngles.size() 
                  << " angles in the system " << std::endl;
        
        std::cout << "These angles are : " << std::endl;
        
        for (std::map<int, std::vector<std::vector<int> > >::iterator tAG 
               =allAnglesIdxs.begin(); tAG != allAnglesIdxs.end(); tAG++)
        {
            std::cout << "There are " << (int)tAG->second.size() 
            << " Angles centered on atom " << allAtoms[tAG->first].id << std::endl
                    << "They are " << std::endl;
                    
            for(std::vector<std::vector<int> >::iterator iAN =tAG->second.begin();
                    iAN != tAG->second.end(); iAN++)
            {
                std::cout << "atoms: " << allAtoms[tAG->first].id <<",  ";
                for (std::vector<int>::iterator iAt = iAN->begin(); 
                        iAt !=iAN->end(); iAt++)
                {
                  std::cout << allAtoms[*iAt].id << ", ";
                }
                std::cout<<std::endl;
            }
        }
    }
    
    void DictCifFile::setOneTorsion(std::vector<int> tAV, REAL tValue,
                                    int tPeriod)
    {
        
        TorsionDict tTor;
      
        tTor.seriNum = itsCurTorsionSeriNum;
        for (int i=0; i <(int)tAV.size(); i++) 
        {
            tTor.atoms.push_back(tAV[i]);
        }
        
        tTor.value  = tValue;
        tTor.period = tPeriod;
        
        allTorsions.push_back(tTor);
        itsCurTorsionSeriNum++;
        
    }
    
    void DictCifFile::SetOneSP2SP2Bond(int tIdx1, int tIdx2)
    {
            // two sp2 atoms
            std::vector<int> tV1, tV2;
            
            int tS1=-1, tS2=-1;
        
            for (std::vector<int>::iterator iAt1=allAtoms[tIdx1].connAtoms.begin();
               iAt1 != allAtoms[tIdx1].connAtoms.end(); iAt1++)
            {
                tS1 =-1;
                for (std::vector<int>::iterator iAt2=allAtoms[tIdx2].connAtoms.begin();
                       iAt2 != allAtoms[tIdx2].connAtoms.end(); iAt2++)
                {
                    tS2 =-1;
                    if (*iAt1 != tIdx2 && *iAt2 !=tIdx1)
                    {
                        if (AtomsInSameRing(allAtoms[*iAt1], allAtoms[*iAt2], allRingsV))
                        {
                            tS1 = *iAt1;
                            tS2 = *iAt2;
                            tV1.push_back(*iAt1);
                            tV2.push_back(*iAt2);
                            std::cout << "atom " << allAtoms[*iAt1].id << " and atom "
                                      << allAtoms[*iAt2].id  << " is in the same ring " 
                                      << std::endl;
                            break;
                        }
                    }
                }
                
                if(tS1 !=-1 && tS2 !=-1)
                {
                    break;
                }
            }
            
            for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
                    iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
            {
                if (*iA1 != tIdx2 && *iA1 !=tS1)
                {
                    tV1.push_back(*iA1);
                }
            }
            
            for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
                    iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
            {
                if (*iA2 != tIdx1 && *iA2 !=tS2)
                {
                    tV2.push_back(*iA2);
                }
            }
            
            if (tS1 ==-1 && tS2 ==-1 
                && (int)tV1.size()==2 && (int)tV2.size()==2)
            {
                int tm;
                if((int)allAtoms[tV1[0]].connAtoms.size() < (int)allAtoms[tV1[1]].connAtoms.size())
                {
                    tm  =tV1[0];
                    tV1[0] = tV1[1];
                    tV1[1] = tm;
                }
                if((int)allAtoms[tV2[0]].connAtoms.size() < (int)allAtoms[tV2[1]].connAtoms.size())
                {
                    tm  =tV2[0];
                    tV2[0] = tV2[1];
                    tV2[1] = tm;
                }
            }
            
            
            REAL va[2][2];
            int per = 2;
            if (tS1 !=-1 && tS2 !=-1)
            {
                va[0][0] =     0.0;
                va[0][1] =   180.0;
                va[1][0] =   180.0;
                va[1][1] =     0.0;
            }
            else
            {
                va[0][0] =  180.0;
                va[0][1] =    0.0;
                va[1][0] =    0.0;
                va[1][1] =  180.0;
            }
            for (int i =0; i < (int)tV1.size(); i++)
            {
                for (int j=0; j < (int)tV2.size(); j++)
                {
                    std::vector<int> aTS;
                    aTS.push_back(tV1[i]);
                    aTS.push_back(tIdx1);
                    aTS.push_back(tIdx2);
                    aTS.push_back(tV2[j]);
                    setOneTorsion(aTS, va[i][j], per);
                }
            } 
            
            /*
            std::cout << "id1 " << allAtoms[tIdx1].id 
                  << " and id2 " << allAtoms[tIdx2].id << std::endl;
            if ((allAtoms[tIdx1].id=="C14" && allAtoms[tIdx2].id=="C23")
                 || (allAtoms[tIdx2].id=="C23" && allAtoms[tIdx1].id=="C14"))
            {
                std::cout << "tS1 = " << tS1 << std::endl
                          << "tS2 = " << tS2 << std::endl;
                
            }
             */
            std::cout << "excluded " << allAtoms[tIdx2].id << ", atom "
                << allAtoms[tIdx1].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV1.begin(); iA !=tV1.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
                
        std::cout << "excluded " << allAtoms[tIdx1].id << " and atom "
                << allAtoms[tIdx2].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
    }
    
    void DictCifFile::SetOneSP2SP3Bond(int tIdx1, int tIdx2, std::string tF)
    {
                
        
        
        // one sp2 atom and one sp3 atom
        std::vector<int> tV1, tV2;
        
        int tS1=-1, tS2=-1;
        
        for (std::vector<int>::iterator iAt1=allAtoms[tIdx1].connAtoms.begin();
               iAt1 != allAtoms[tIdx1].connAtoms.end(); iAt1++)
        {
            tS1 =-1;
            for (std::vector<int>::iterator iAt2=allAtoms[tIdx2].connAtoms.begin();
                       iAt2 != allAtoms[tIdx2].connAtoms.end(); iAt2++)
            {
                tS2 =-1;
                if (*iAt1 != tIdx2 && *iAt2 !=tIdx1)
                {
                    if (AtomsInSameRing(allAtoms[*iAt1], allAtoms[*iAt2], allRingsV))
                    {
                            tS1 = *iAt1;
                            tS2 = *iAt2;
                            tV1.push_back(*iAt1);
                            tV2.push_back(*iAt2);
                            std::cout << "atom " << allAtoms[*iAt1].id << " and atom "
                                      << allAtoms[*iAt2].id  << " is in the same ring " 
                                      << std::endl;
                            break;
                            
                    }
                }
            }
            if(tS1 !=-1 && tS2 !=-1)
            {
                
                break;
                
            }
        }
            
        REAL va1  =   0.0;
        REAL va2  =  60.0;
        int  per  =     6;
        
        REAL va[2][3];
        if (tS1 !=-1 && tS2 !=-1)
        {
            va[0][0] =    va1;
            va[0][1] =  2*va2;
            va[0][2] = -2*va2;
        
            va[1][0] =  3*va2;
            va[1][1] =  -va2;
            va[1][2] =   va2;
            /*
            if (tF=="even")
            {
                va[0][0] =    va2;
                va[0][1] =  3*va2;
                va[0][2] =   -va2;
        
                va[1][0] =  -2*va2;
                va[1][1] =     va1;
                va[1][2] =   2*va2;
            }
            else
            {
                va[0][0] =   -va2;
                va[0][1] =    va2;
                va[0][2] =   3*va2;
        
                va[1][0] =   2*va2;
                va[1][1] =  -2*va2;
                va[1][2] =     va1;
            }
             */
        }
        else
        {  /*
            va[0][0] =    va1;
            va[0][1] =  2*va2;
            va[0][2] = -2*va2;
        
            va[1][0] =  3*va2;
            va[1][1] =  -va2;
            va[1][2] =   va2;
           */  
            
            va[0][0] =  3*va2;
            va[0][1] =  -va2;
            va[0][2] =   va2;
            
            va[1][0] =    va1;
            va[1][1] =  2*va2;
            va[1][2] = -2*va2;
        }
        
        
        for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
        {
            if (*iA1 != tIdx2 && *iA1 != tS1)
            {
                tV1.push_back(*iA1);
            }
        }
        
        // for sp3 atoms
        if ((int)allAtoms[tIdx2].inChirals.size() !=0)
        {
            
            int iCh = allAtoms[tIdx2].inChirals[0];
            
            
            /*
            if(tS1 !=-1 && tS2 !=-1)
            {
                std::cout << "Chiral center " << allAtoms[tIdx2].id << std::endl;
                std::cout << "Chiral atoms\n";
                for (int i=0; i < (int)allChirals[iCh].atoms.size(); i++)
                {
                    std::cout << allAtoms[allChirals[iCh].atoms[i]].id << std::endl;
                }
                exit(1);
                        
            }
            */
            // buildChiralCluster2(allChirals[iCh], tV2, tIdx1, allAtoms[tIdx2].connAtoms);
            buildChiralCluster2(allChirals[iCh], tV2, tIdx1);
           /*
            //if(allAtoms[tIdx2].id =="C15" || allAtoms[tIdx2].id =="C9")
            //{
             std::cout << "Chiral serial number is " << iCh << std::endl;   
             std::cout << "root atom is  " << allAtoms[tIdx1].id << std::endl;
             std::cout << "atom  " << allAtoms[tIdx2].id << " mutable size is "
                      << (int)allChirals[iCh].mutTable.size() << std::endl
                      << std::endl << allChirals[iCh].mutTable[tIdx1][0]
                      << " and " << allChirals[iCh].mutTable[tIdx1][1] << std::endl
                      << " chiral atom seq is " << std::endl; 
             for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
             {
                 std::cout << "atom " << allAtoms[*iA].id << std::endl;
             }
           
            //}
            */
        }
        else
        {
            for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
            {
                if (*iA2 !=tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end() &&
                        allAtoms[*iA2].chemType=="H")
                {
                    tV2.push_back(*iA2);
                    break;
                }
            }
        }
        
        // put those atoms which are not, such as H atoms, in chiral 
        // list into the output list
        
        for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
        {
            if (*iA2 != tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end())
            {
                tV2.push_back(*iA2);
            }
        }
        
        if ((int)tV1.size() <=(int)tV2.size())
        {
            for (int i =0; i < (int)tV1.size(); i++)
            {
                for (int j=0; j < (int)tV2.size(); j++)
                {
                    std::vector<int> aTS;
                    aTS.push_back(tV1[i]);
                    aTS.push_back(tIdx1);
                    aTS.push_back(tIdx2);
                    aTS.push_back(tV2[j]);
                    setOneTorsion(aTS, va[i][j], per);
                }
            }
        }
        else
        {
            for (int i =0; i < (int)tV2.size(); i++)
            {
                for (int j=0; j < (int)tV1.size(); j++)
                {
                    std::vector<int> aTS;
                    aTS.push_back(tV2[i]);
                    aTS.push_back(tIdx1);
                    aTS.push_back(tIdx2);
                    aTS.push_back(tV1[j]);
                    setOneTorsion(aTS, va[i][j], per);
                }
            }   
        }
        
        std::cout << "excluded " << allAtoms[tIdx2].id << ", atom "
                << allAtoms[tIdx1].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV1.begin(); iA !=tV1.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
                
        std::cout << "excluded " << allAtoms[tIdx1].id << " and atom "
                << allAtoms[tIdx2].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
        
       
    }
    
    
    void DictCifFile::SetOneSP2SP3Bond(int tIdx1, int tIdx2)
    {
                
        
        
        // one sp2 atom and one sp3 atom
        std::vector<int> tV1, tV2;
        
        int tS1=-1, tS2=-1;
        
        for (std::vector<int>::iterator iAt1=allAtoms[tIdx1].connAtoms.begin();
               iAt1 != allAtoms[tIdx1].connAtoms.end(); iAt1++)
        {
            tS1 =-1;
            for (std::vector<int>::iterator iAt2=allAtoms[tIdx2].connAtoms.begin();
                       iAt2 != allAtoms[tIdx2].connAtoms.end(); iAt2++)
            {
                tS2 =-1;
                if (*iAt1 != tIdx2 && *iAt2 !=tIdx1)
                {
                    if (AtomsInSameRing(allAtoms[*iAt1], allAtoms[*iAt2], allRingsV))
                    {
                            tS1 = *iAt1;
                            tS2 = *iAt2;
                            tV1.push_back(*iAt1);
                            tV2.push_back(*iAt2);
                            std::cout << "atom " << allAtoms[*iAt1].id << " and atom "
                                      << allAtoms[*iAt2].id  << " is in the same ring " 
                                      << std::endl;
                            break;
                            
                    }
                }
            }
            if(tS1 !=-1 && tS2 !=-1)
            {
                
                break;
                
            }
        }
            
        REAL va1  =   0.0;
        REAL va2  =  60.0;
        int  per  =     6;
        
        REAL va[2][3];
        if (tS1 !=-1 && tS2 !=-1)
        {
            va[0][0] =    va1;
            va[0][1] =  2*va2;
            va[0][2] = -2*va2;
        
            va[1][0] =  3*va2;
            va[1][1] =  -va2;
            va[1][2] =   va2;
            /*
            va[0][0] =    va2;
            va[0][1] =  3*va2;
            va[0][2] =   -va2;
        
            va[1][0] =  -2*va2;
            va[1][1] =     va1;
            va[1][2] =   2*va2;
             */
        }
        else
        {
           /*
            va[0][0] =    va1;
            va[0][1] =  2*va2;
            va[0][2] = -2*va2;
        
            va[1][0] =  3*va2;
            va[1][1] =  -va2;
            va[1][2] =   va2;
           */
           
            va[0][0] =  3*va2;
            va[0][1] =  -va2;
            va[0][2] =   va2;
            
            va[1][0] =    va1;
            va[1][1] =  2*va2;
            va[1][2] = -2*va2;
        }
        
        
        for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
        {
            if (*iA1 != tIdx2 && *iA1 != tS1)
            {
                tV1.push_back(*iA1);
            }
        }
        
        // for sp3 atoms
        if ((int)allAtoms[tIdx2].inChirals.size() !=0)
        {
            
            int iCh = allAtoms[tIdx2].inChirals[0];
            
            
            /*
            if(tS1 !=-1 && tS2 !=-1)
            {
                std::cout << "Chiral center " << allAtoms[tIdx2].id << std::endl;
                std::cout << "Chiral atoms\n";
                for (int i=0; i < (int)allChirals[iCh].atoms.size(); i++)
                {
                    std::cout << allAtoms[allChirals[iCh].atoms[i]].id << std::endl;
                }
                exit(1);
                        
            }
            */
            // buildChiralCluster2(allChirals[iCh], tV2, tIdx1, allAtoms[tIdx2].connAtoms);
            buildChiralCluster2(allChirals[iCh], tV2, tIdx1);
           /*
            //if(allAtoms[tIdx2].id =="C15" || allAtoms[tIdx2].id =="C9")
            //{
             std::cout << "Chiral serial number is " << iCh << std::endl;   
             std::cout << "root atom is  " << allAtoms[tIdx1].id << std::endl;
             std::cout << "atom  " << allAtoms[tIdx2].id << " mutable size is "
                      << (int)allChirals[iCh].mutTable.size() << std::endl
                      << std::endl << allChirals[iCh].mutTable[tIdx1][0]
                      << " and " << allChirals[iCh].mutTable[tIdx1][1] << std::endl
                      << " chiral atom seq is " << std::endl; 
             for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
             {
                 std::cout << "atom " << allAtoms[*iA].id << std::endl;
             }
           
            //}
            */
        }
        else
        {
            for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
            {
                if (*iA2 !=tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end() &&
                        allAtoms[*iA2].chemType=="H")
                {
                    tV2.push_back(*iA2);
                    break;
                }
            }
        }
        
        // put those atoms which are not, such as H atoms, in chiral 
        // list into the output list
        
        for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
        {
            if (*iA2 != tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end())
            {
                tV2.push_back(*iA2);
            }
        }
        
        if ((int)tV1.size() <=(int)tV2.size())
        {
            for (int i =0; i < (int)tV1.size(); i++)
            {
                for (int j=0; j < (int)tV2.size(); j++)
                {
                    std::vector<int> aTS;
                    aTS.push_back(tV1[i]);
                    aTS.push_back(tIdx1);
                    aTS.push_back(tIdx2);
                    aTS.push_back(tV2[j]);
                    setOneTorsion(aTS, va[i][j], per);
                }
            }
        }
        else
        {
            for (int i =0; i < (int)tV2.size(); i++)
            {
                for (int j=0; j < (int)tV1.size(); j++)
                {
                    std::vector<int> aTS;
                    aTS.push_back(tV2[i]);
                    aTS.push_back(tIdx1);
                    aTS.push_back(tIdx2);
                    aTS.push_back(tV1[j]);
                    setOneTorsion(aTS, va[i][j], per);
                }
            }   
        }
        
        std::cout << "excluded " << allAtoms[tIdx2].id << ", atom "
                << allAtoms[tIdx1].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV1.begin(); iA !=tV1.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
                
        std::cout << "excluded " << allAtoms[tIdx1].id << " and atom "
                << allAtoms[tIdx2].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
        
       
    }
    
    
    
    void DictCifFile::SetOneSP3SP3Bond(int tIdx1, int tIdx2)
    {
        
        REAL va[3][3];
        va[0][0] = 60.0;
        va[0][1] = 180.0;
        va[0][2] = -60.0;
        
        va[1][0] = -60.0;
        va[1][1] =  60.0;
        va[1][2] = 180.0;
        
        va[2][0] = 180.0;
        va[2][1] = -60.0;
        va[2][2] =  60.0;
        
        int  per =     3;
        
        // One sp3 atom and one sp3 atom
        // Several procedures
        // 1. if there are two atoms in the same ring
        
        std::vector<int> tV1, tV2;
        
        int tS1=-1, tS2=-1;
        
        for (std::vector<int>::iterator iAt1=allAtoms[tIdx1].connAtoms.begin();
               iAt1 != allAtoms[tIdx1].connAtoms.end(); iAt1++)
        {
            tS1 =-1;
            for (std::vector<int>::iterator iAt2=allAtoms[tIdx2].connAtoms.begin();
                       iAt2 != allAtoms[tIdx2].connAtoms.end(); iAt2++)
            {
                tS2 =-1;
                if (*iAt1 != tIdx2 && *iAt2 !=tIdx1)
                {
                    if (AtomsInSameRing(allAtoms[*iAt1], allAtoms[*iAt2], allRingsV))
                    {
                            tS1 = *iAt1;
                            tS2 = *iAt2;
                            tV1.push_back(*iAt1);
                            tV2.push_back(*iAt2);
                            std::cout << "atom " << allAtoms[*iAt1].id << " and atom "
                                      << allAtoms[*iAt2].id  << " is in the same ring " 
                                      << std::endl;
                            break;
                    }
                }
            }
            if(tS1 !=-1 && tS2 !=-1)
            {
                break;
            }
        }
        
        // 2. take consider of chirals 
        // std::cout << "tS1 " << tS1 << " tS2 " << tS2 << std::endl;
        
        if ((int)allAtoms[tIdx1].inChirals.size() !=0)
        {
            int iCh = allAtoms[tIdx1].inChirals[0];
            
            // buildChiralCluster2(allChirals[iCh], tV1, tIdx2, allAtoms[tIdx1].connAtoms);
            buildChiralCluster2(allChirals[iCh], tV1, tIdx2);
        }
        else
        {
            for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
            {
                if (*iA1 !=tIdx2 && std::find(tV1.begin(), tV1.end(), *iA1)==tV1.end() &&
                        allAtoms[*iA1].chemType=="H")
                {
                    tV1.push_back(*iA1);
                    break;
                }
            }
        }
        
       
        if ((int)allAtoms[tIdx2].inChirals.size() !=0)
        {
            int iCh = allAtoms[tIdx2].inChirals[0];
            
            //buildChiralCluster2(allChirals[iCh], tV2, tIdx1, allAtoms[tIdx2].connAtoms);
            buildChiralCluster2(allChirals[iCh], tV2, tIdx1);
        }
        else
        {
            for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
            {
                if (*iA2 !=tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end() &&
                        allAtoms[*iA2].chemType=="H")
                {
                    tV2.push_back(*iA2);
                    break;
                }
            }
        }
       
        
       
           
        // the rest atoms not included in the chiral atom cluster 
        for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
        {
            if (*iA1 !=tIdx2 && std::find(tV1.begin(), tV1.end(), *iA1)==tV1.end())
            {
                tV1.push_back(*iA1);
            }
        }
        
        for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
        {
            if (*iA2 !=tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end())
            {
                tV2.push_back(*iA2);
            }
        }
        
        std::cout << "excluded " << allAtoms[tIdx2].id << " and atom "
                << allAtoms[tIdx1].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV1.begin(); iA !=tV1.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
                
        std::cout << "excluded " << allAtoms[tIdx1].id << " and atom "
                << allAtoms[tIdx2].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
        
        
        for (int i =0; i < (int)tV1.size(); i++)
        {
            for (int j=0; j < (int)tV2.size(); j++)
            {
                std::vector<int> aTS;
                aTS.push_back(tV1[i]);
                aTS.push_back(tIdx1);
                aTS.push_back(tIdx2);
                aTS.push_back(tV2[j]);
                setOneTorsion(aTS, va[i][j], per);
            }
        }
    }
    
    void DictCifFile::SetOneSP3SP3Bond(int tIdx1, int tIdx2, std::string tF)
    {
        
        REAL va[3][3];
        if (tF=="even")
        {
            va[0][0] =  60.0;
            va[0][1] = 180.0;
            va[0][2] = -60.0;
        
            va[1][0] = -60.0;
            va[1][1] =  60.0;
            va[1][2] = 180.0;
        
            va[2][0] = 180.0;
            va[2][1] = -60.0;
            va[2][2] =  60.0;
        }
        else if(tF=="odd")
        {
            va[0][0] = -60.0;
            va[0][1] =  60.0;
            va[0][2] = 180.0;
        
            va[1][0] = 180.0;
            va[1][1] = -60.0;
            va[1][2] =  60.0;
        
            va[2][0] =  60.0;
            va[2][1] = 180.0;
            va[2][2] = -60.0;
        }
        else
        {
            std::cout << "what is the sequence idx of this torsion in a ring, even or odd "
                    << std::endl;
            exit(1);
        }
        
        int  per =     3;
        
        // One sp3 atom and one sp3 atom
        // Several procedures
        // 1. if there are two atoms in the same ring
        
        std::vector<int> tV1, tV2;
        
        int tS1=-1, tS2=-1;
        
        for (std::vector<int>::iterator iAt1=allAtoms[tIdx1].connAtoms.begin();
               iAt1 != allAtoms[tIdx1].connAtoms.end(); iAt1++)
        {
            tS1 =-1;
            for (std::vector<int>::iterator iAt2=allAtoms[tIdx2].connAtoms.begin();
                       iAt2 != allAtoms[tIdx2].connAtoms.end(); iAt2++)
            {
                tS2 =-1;
                if (*iAt1 != tIdx2 && *iAt2 !=tIdx1)
                {
                    if (AtomsInSameRing(allAtoms[*iAt1], allAtoms[*iAt2], allRingsV))
                    {
                            tS1 = *iAt1;
                            tS2 = *iAt2;
                            tV1.push_back(*iAt1);
                            tV2.push_back(*iAt2);
                            std::cout << "atom " << allAtoms[*iAt1].id << " and atom "
                                      << allAtoms[*iAt2].id  << " is in the same ring " 
                                      << std::endl;
                            break;
                    }
                }
            }
            if(tS1 !=-1 && tS2 !=-1)
            {
                break;
            }
        }
        
        // 2. take consider of chirals 
        // std::cout << "tS1 " << tS1 << " tS2 " << tS2 << std::endl;
        
        if ((int)allAtoms[tIdx1].inChirals.size() !=0)
        {
            int iCh = allAtoms[tIdx1].inChirals[0];
            
            // buildChiralCluster2(allChirals[iCh], tV1, tIdx2, allAtoms[tIdx1].connAtoms);
            buildChiralCluster2(allChirals[iCh], tV1, tIdx2);
        }
        else
        {
            for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
            {
                if (*iA1 !=tIdx2 && std::find(tV1.begin(), tV1.end(), *iA1)==tV1.end() &&
                        allAtoms[*iA1].chemType=="H")
                {
                    tV1.push_back(*iA1);
                    break;
                }
            }
        }
        
       
        if ((int)allAtoms[tIdx2].inChirals.size() !=0)
        {
            int iCh = allAtoms[tIdx2].inChirals[0];
            
            //buildChiralCluster2(allChirals[iCh], tV2, tIdx1, allAtoms[tIdx2].connAtoms);
            buildChiralCluster2(allChirals[iCh], tV2, tIdx1);
        }
        else
        {
            for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
            {
                if (*iA2 !=tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end() &&
                        allAtoms[*iA2].chemType=="H")
                {
                    tV2.push_back(*iA2);
                    break;
                }
            }
        }
        
           
        // the rest atoms not included in the chiral atom cluster 
        for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
        {
            if (*iA1 !=tIdx2 && std::find(tV1.begin(), tV1.end(), *iA1)==tV1.end())
            {
                tV1.push_back(*iA1);
            }
        }
        
        for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
        {
            if (*iA2 !=tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end())
            {
                tV2.push_back(*iA2);
            }
        }
        
        std::cout << "excluded " << allAtoms[tIdx2].id << " and atom "
                << allAtoms[tIdx1].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV1.begin(); iA !=tV1.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
                
        std::cout << "excluded " << allAtoms[tIdx1].id << " and atom "
                << allAtoms[tIdx2].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
        
        
        for (int i =0; i < (int)tV1.size(); i++)
        {
            for (int j=0; j < (int)tV2.size(); j++)
            {
                std::vector<int> aTS;
                aTS.push_back(tV1[i]);
                aTS.push_back(tIdx1);
                aTS.push_back(tIdx2);
                aTS.push_back(tV2[j]);
                setOneTorsion(aTS, va[i][j], per);
            }
        }
    }
    
    void DictCifFile::SetOneSP3SP3Bond4H(int tIdx1, int tIdx2)
    {
        // Each sp3 center atoms connected two H atoms
        // Two sp3 center atoms are not in the same ring
        
        REAL va[3][3];
        va[0][0] = 180.0;
        va[0][1] = -60.0;
        va[0][2] =  60.0;
        
        va[1][0] =  60.0;
        va[1][1] =  180.0;
        va[1][2] = -60.0;
        
        va[2][0] = -60.0;
        va[2][1] =  60.0;
        va[2][2] = 180.0;
        
        int  per =     3;
        
        std::vector<int> tV1, tV2;
        
        for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
        {
            if (*iA1 !=tIdx2 && std::find(tV1.begin(), tV1.end(), *iA1)==tV1.end() &&
                        allAtoms[*iA1].chemType !="H")
            {
                tV1.push_back(*iA1);
                break;
            }
        }
        
        for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
        {
            if (*iA2 !=tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end() &&
                        allAtoms[*iA2].chemType !="H")
            {
                tV2.push_back(*iA2);
                break;
            }
        }
        
        // the rest atoms not included in the chiral atom cluster 
        for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
        {
            if (*iA1 !=tIdx2 && std::find(tV1.begin(), tV1.end(), *iA1)==tV1.end())
            {
                tV1.push_back(*iA1);
            }
        }
        
        for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
        {
            if (*iA2 !=tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end())
            {
                tV2.push_back(*iA2);
            }
        }
        
        std::cout << "excluded " << allAtoms[tIdx2].id << " and atom "
                << allAtoms[tIdx1].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV1.begin(); iA !=tV1.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
                
        std::cout << "excluded " << allAtoms[tIdx1].id << " and atom "
                << allAtoms[tIdx2].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
        
        
        for (int i =0; i < (int)tV1.size(); i++)
        {
            for (int j=0; j < (int)tV2.size(); j++)
            {
                std::vector<int> aTS;
                aTS.push_back(tV1[i]);
                aTS.push_back(tIdx1);
                aTS.push_back(tIdx2);
                aTS.push_back(tV2[j]);
                setOneTorsion(aTS, va[i][j], per);
            }
        }
   
    }
    
    bool DictCifFile::checkSP3SP34H(int tIdx1, int tIdx2)
    {
        int  n1  = 0;
        int  n2  = 0;
        bool h4  = false;
        
        for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
        {
            if (*iA1 !=tIdx2 && allAtoms[*iA1].chemType =="H")
            {
                n1++;
            }
        }
        
        for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
               iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
        {
            if (*iA2!=tIdx1 && allAtoms[*iA2].chemType =="H")
            {
                n2++;
            }
        }
        
        if (n1==2 && n2==2)
        {
            h4 = true;
        }
        
        return h4;
    }
    
    void DictCifFile::SetOneSP3OxyColumnBond(int tIdx1, int tIdx2, 
                                             int tPer, REAL tIniValue)
    {
        // atoms connected to atoms tIdx1, tIdx2
        int iS=-1;
        std::vector<int> tV1, tV2;
        int tS1=-1, tS2=-1;
        // std::cout << "For atom " << allAtoms[tIdx1].id << " and atom " 
        //          << allAtoms[tIdx2].id << "connected atoms are " << std::endl;
        
        for (std::vector<int>::iterator iAt1=allAtoms[tIdx1].connAtoms.begin();
               iAt1 != allAtoms[tIdx1].connAtoms.end(); iAt1++)
        {
            tS1 =-1;
            // std::cout << "1 linked " << allAtoms[*iAt1].id << std::endl;
            
            for (std::vector<int>::iterator iAt2=allAtoms[tIdx2].connAtoms.begin();
                       iAt2 != allAtoms[tIdx2].connAtoms.end(); iAt2++)
            {
                tS2 =-1;
                //std::cout << "2 linked " << allAtoms[*iAt2].id << std::endl;
                if (*iAt1 != tIdx2 && *iAt2 !=tIdx1)
                {
                    if (AtomsInSameRing(allAtoms[*iAt1], allAtoms[*iAt2], allRingsV))
                    {
                            tS1 = *iAt1;
                            tS2 = *iAt2;
                            tV1.push_back(*iAt1);
                            tV2.push_back(*iAt2);
                            //std::cout << "atom " << allAtoms[*iAt1].id << " and atom "
                            //          << allAtoms[*iAt2].id  << " is in the same ring " 
                            //          << std::endl;
                            break;
                    }
                }
            }
            if(tS1 !=-1 && tS2 !=-1)
            {
                break;
            }
        }
        
        if(tS1 ==-1 && tS2 ==-1)
        {
            //First find a atom which is not H as a start atom
            
            for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
                    iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
            {
                if(allAtoms[*iA1].chemType !="H" && *iA1 !=tIdx2)
                {
                    iS=*iA1;
                    break;
                }
            }
            
            // if in some cases, only H has involved, use them.
            if (iS==-1)
            {
                iS=allAtoms[tIdx1].connAtoms[0];
            }
            
            tV1.push_back(iS);
        }
        
            if ((int)allAtoms[tIdx1].inChirals.size() !=0)
            {
                int iCh = allAtoms[tIdx1].inChirals[0];
                //buildChiralCluster2(allChirals[iCh], tV1, tIdx2, allAtoms[tIdx1].connAtoms);
                //std::cout << "buildChiralCluster2 " << std::endl;
                buildChiralCluster2(allChirals[iCh], tV1, tIdx2);
            }
            
        //std::cout << "tV1 " << (int)tV1.size() << std::endl;
        for(std::vector<int>::iterator iV1=tV1.begin(); iV1 != tV1.end();
                iV1++)
        {
            std::cout << allAtoms[*iV1].id  << std::endl;
        }
           
                
            for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
            {
                if (*iA1 !=tIdx2 && std::find(tV1.begin(), tV1.end(), *iA1)==tV1.end())
                {
                    tV1.push_back(*iA1);
                }
            }
           
            for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
                    iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
            {
                if (*iA2 != tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end())
                {
                    tV2.push_back(*iA2);
                }
            }
            
           // std::cout << "number of tV2 " << (int)tV2.size() << std::endl;
            
            REAL perValue = 360.0/tPer;
            REAL iniValue;
            
            int n = 0;
            
            for (int i=0; i <(int)tV1.size(); i++)
            {
                if(tS1 !=-1 && tS2 !=-1)
                {
                   iniValue= 60.0 +  n*perValue;
                }
                else
                {
                    iniValue=tIniValue + n*perValue;
                }
                // std::cout << "iniValue " << iniValue << std::endl;
                
                int m =0;
                for (int j=0; j < (int)tV2.size(); j++)
                {
                    //std::cout << "tV1 " << i << " value " << allAtoms[tV1[i]].id << std::endl;
                    //std::cout << "tV2 " << j << " value " << allAtoms[tV2[j]].id << std::endl;
                    std::vector<int> TS;
                    TS.push_back(tV1[i]);
                    TS.push_back(tIdx1);
                    TS.push_back(tIdx2);
                    TS.push_back(tV2[j]);
                    
                 
                    REAL curValue = iniValue +m*perValue;
                  
                    if (curValue > 360.0)
                    {
                        curValue = curValue -360.0;
                    }
                    else if (curValue > 180.0)
                    {
                        curValue = curValue -360.0;
                    }
                    setOneTorsion(TS, curValue, tPer);
                    //std::cout << "atom 1 " << allAtoms[tV1[i]].id
                    //          << " atom 4 " << allAtoms[tV2[j]].id
                    //          << " tor " << curValue << std::endl;
                    m++;
                }
                n++;
            }
            std::cout << "excluded " << allAtoms[tIdx2].id << ", atom "
                << allAtoms[tIdx1].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV1.begin(); iA !=tV1.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
                
        std::cout << "excluded " << allAtoms[tIdx1].id << " and atom "
                << allAtoms[tIdx2].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
    }
    
    void DictCifFile::SetOneSP3OxyColumnBond(int tIdx1, int tIdx2, 
                                             int tPer, REAL tIniValue, 
                                             std::string tF)
    {
        // atoms connected to atoms tIdx1, tIdx2
        int iS=-1;
        std::vector<int> tV1, tV2;
        int tS1=-1, tS2=-1;
        //std::cout << "For atom " << allAtoms[tIdx1].id << " and atom " 
        //          << allAtoms[tIdx2].id << "connected atoms are " << std::endl;
        
        for (std::vector<int>::iterator iAt1=allAtoms[tIdx1].connAtoms.begin();
               iAt1 != allAtoms[tIdx1].connAtoms.end(); iAt1++)
        {
            tS1 =-1;
            //std::cout << "1 linked " << allAtoms[*iAt1].id << std::endl;
            
            for (std::vector<int>::iterator iAt2=allAtoms[tIdx2].connAtoms.begin();
                       iAt2 != allAtoms[tIdx2].connAtoms.end(); iAt2++)
            {
                tS2 =-1;
                //std::cout << "2 linked " << allAtoms[*iAt2].id << std::endl;
                if (*iAt1 != tIdx2 && *iAt2 !=tIdx1)
                {
                    if (AtomsInSameRing(allAtoms[*iAt1], allAtoms[*iAt2], allRingsV))
                    {
                            tS1 = *iAt1;
                            tS2 = *iAt2;
                            tV1.push_back(*iAt1);
                            tV2.push_back(*iAt2);
                            std::cout << "atom " << allAtoms[*iAt1].id << " and atom "
                                      << allAtoms[*iAt2].id  << " is in the same ring " 
                                      << std::endl;
                            break;
                    }
                }
            }
            if(tS1 !=-1 && tS2 !=-1)
            {
                break;
            }
        }
        
        if(tS1 ==-1 && tS2 ==-1)
        {
            //First find a atom which is not H as a start atom
            
            for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
                    iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
            {
                if(allAtoms[*iA1].chemType !="H" && *iA1 !=tIdx2)
                {
                    iS=*iA1;
                    break;
                }
            }
            
            // if in some cases, only H has involved, use them.
            if (iS==-1)
            {
                iS=allAtoms[tIdx1].connAtoms[0];
            }
            
            tV1.push_back(iS);
        }
        
            if ((int)allAtoms[tIdx1].inChirals.size() !=0)
            {
                int iCh = allAtoms[tIdx1].inChirals[0];
                //buildChiralCluster2(allChirals[iCh], tV1, tIdx2, allAtoms[tIdx1].connAtoms);
                //std::cout << "buildChiralCluster2 " << std::endl;
                buildChiralCluster2(allChirals[iCh], tV1, tIdx2);
            }
            
        // std::cout << "tV1 " << (int)tV1.size() << std::endl;
        for(std::vector<int>::iterator iV1=tV1.begin(); iV1 != tV1.end();
                iV1++)
        {
            std::cout << allAtoms[*iV1].id  << std::endl;
        }
            
            for (std::vector<int>::iterator iA1=allAtoms[tIdx1].connAtoms.begin();
               iA1 != allAtoms[tIdx1].connAtoms.end(); iA1++)
            {
                if (*iA1 !=tIdx2 && std::find(tV1.begin(), tV1.end(), *iA1)==tV1.end())
                {
                    tV1.push_back(*iA1);
                }
            }
           
            for (std::vector<int>::iterator iA2=allAtoms[tIdx2].connAtoms.begin();
                    iA2 != allAtoms[tIdx2].connAtoms.end(); iA2++)
            {
                if (*iA2 != tIdx1 && std::find(tV2.begin(), tV2.end(), *iA2)==tV2.end())
                {
                    tV2.push_back(*iA2);
                }
            }
            
           // std::cout << "number of tV2 " << (int)tV2.size() << std::endl;
            
            REAL perValue = 360.0/tPer;
            REAL iniValue = 0.0;
            
            int n = 0;
            
            for (int i=0; i <(int)tV1.size(); i++)
            {
                if(tS1 !=-1 && tS2 !=-1)
                {
                    if (tF=="even")
                    {
                        iniValue= 60.0 +  n*perValue;
                    }
                    else if (tF=="odd")
                    {
                        iniValue= -60.0 +  n*perValue;
                    }
                }
                else
                {
                    iniValue=tIniValue + n*perValue;
                }
                std::cout << "iniValue " << iniValue << std::endl;
                
                int m =0;
                for (int j=0; j < (int)tV2.size(); j++)
                {
                    std::cout << "tV1 " << i << " value " << allAtoms[tV1[i]].id << std::endl;
                    std::cout << "tV2 " << j << " value " << allAtoms[tV2[j]].id << std::endl;
                    std::vector<int> TS;
                    TS.push_back(tV1[i]);
                    TS.push_back(tIdx1);
                    TS.push_back(tIdx2);
                    TS.push_back(tV2[j]);
                    
                 
                    REAL curValue = iniValue +m*perValue;
                  
                    if (curValue > 360.0)
                    {
                        curValue = curValue -360.0;
                    }
                    else if (curValue > 180.0)
                    {
                        curValue = curValue -360.0;
                    }
                    setOneTorsion(TS, curValue, tPer);
                    std::cout << "atom 1 " << allAtoms[tV1[i]].id
                              << " atom 4 " << allAtoms[tV2[j]].id
                              << " tor " << curValue << std::endl;
                    m++;
                }
                n++;
            }
            std::cout << "excluded " << allAtoms[tIdx2].id << ", atom "
                << allAtoms[tIdx1].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV1.begin(); iA !=tV1.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
                
        std::cout << "excluded " << allAtoms[tIdx1].id << " and atom "
                << allAtoms[tIdx2].id << " form torsion using " << std::endl;
        for (std::vector<int>::iterator iA=tV2.begin(); iA !=tV2.end(); iA++)
        {
            std::cout << "atom " << allAtoms[*iA].id << std::endl;
        }
    }
    
    
    void DictCifFile::setTorsionFromOneBond(int tIdx1, int tIdx2)
    {
        int bIdx1 =  allAtoms[tIdx1].bondingIdx;
        ID  id1   =  allAtoms[tIdx1].chemType;
        
        int bIdx2 =  allAtoms[tIdx2].bondingIdx;
        ID  id2   =  allAtoms[tIdx2].chemType;
        
        //std::cout << "id1 " << allAtoms[tIdx1].id 
        //          << " and id2 " << allAtoms[tIdx2].id << std::endl;
        
        std::vector<ID> OxyCol;
        
        OxyCol.push_back("O");
        OxyCol.push_back("S");
        OxyCol.push_back("Se");
        OxyCol.push_back("Te");
        OxyCol.push_back("Po");
        
        std::vector<ID>::iterator iFind1 = std::find(OxyCol.begin(), OxyCol.end(), id1);
        std::vector<ID>::iterator iFind2 = std::find(OxyCol.begin(), OxyCol.end(), id2);
        
        if ((iFind1 != OxyCol.end() && bIdx1==2) && 
                (iFind2 != OxyCol.end() && bIdx2==3))
        {
            // two Oxy column atoms with sp3 orbiting
            //std::cout << "SetOneSP3OxyColumnBond(bIdx1, bIdx2, 2, 90.0)" << std::endl;
            //std::cout << "atom 1 "  << allAtoms[tIdx1].id 
            //          << " atom 2 " << allAtoms[tIdx2].id << std::endl;
            SetOneSP3OxyColumnBond(tIdx2, tIdx1, 3, 90.0);   
        }
        else if ((iFind1 != OxyCol.end() && bIdx1==2) && bIdx2==3)  
        {
            //std::cout << "SetOneSP3OxyColumnBond(bIdx1, bIdx2, 12, 180.0)" << std::endl;
            //std::cout << "atom 1 sp2 " << allAtoms[tIdx1].id 
            //        << " atom 2 sp3 " << allAtoms[tIdx2].id << std::endl;
            SetOneSP3OxyColumnBond(tIdx2, tIdx1, 3, 180.0);
            
        }
        else if ((iFind2 != OxyCol.end() && bIdx2==2) && bIdx1==3)  
        {
            //std::cout << "SetOneSP3OxyColumnBond" << std::endl;
            //std::cout << "atom 1 sp3 " << allAtoms[tIdx1].id 
            //          << "atom 2 sp2 " << allAtoms[tIdx2].id << std::endl;
            SetOneSP3OxyColumnBond(tIdx1, tIdx2, 3, 180.0);
        }
        else if (bIdx1==2 && bIdx2==2)
        {
            //std::cout << "SetOneSP2SP2Bond" << std::endl;
            SetOneSP2SP2Bond(tIdx1, tIdx2);
            
        }
        else if ((bIdx1==2 && bIdx2==3) || (bIdx1==3 && bIdx2==2))
        {
            if (bIdx1==2 && bIdx2==3)
            {
                //std::cout << "SetOneSP2SP3Bond(bIdx1, bIdx2)" << std::endl;
                //std::cout << " atom 1 sp2 " << allAtoms[tIdx1].id 
                //          << " atom 2 sp3 " << allAtoms[tIdx2].id << std::endl;
                SetOneSP2SP3Bond(tIdx1, tIdx2);
            }
            else
            {  
                //std::cout << "SetOneSP2SP3Bond(bIdx1, bIdx2)" << std::endl;
                //std::cout << " atom 1 sp3 " << allAtoms[tIdx1].id 
                //          << " atom 2 sp2 " << allAtoms[tIdx2].id << std::endl;
                SetOneSP2SP3Bond(tIdx2, tIdx1);
            }
        }
        else if (bIdx1==3 && bIdx2==3)
        {
            // std::cout << "SetOneSP3SP3Bond" << std::endl;
            if (!AtomsInSameRing(allAtoms[tIdx1], allAtoms[tIdx2], allRingsV)
                 && checkSP3SP34H(tIdx1, tIdx2))
            {
                SetOneSP3SP3Bond4H(tIdx1, tIdx2);
            }
            else
            {
                SetOneSP3SP3Bond(tIdx1, tIdx2);
            }
        }  
    }
     
    void DictCifFile::setTorsionFromOneBond(int tIdx1, int tIdx2, std::string tF)
    {
                int bIdx1 =  allAtoms[tIdx1].bondingIdx;
        ID  id1   =  allAtoms[tIdx1].chemType;
        
        int bIdx2 =  allAtoms[tIdx2].bondingIdx;
        ID  id2   =  allAtoms[tIdx2].chemType;
        
        //std::cout << "id1 " << allAtoms[tIdx1].id 
        //          << " and id2 " << allAtoms[tIdx2].id << std::endl;
        
        std::vector<ID> OxyCol;
        
        OxyCol.push_back("O");
        OxyCol.push_back("S");
        OxyCol.push_back("Se");
        OxyCol.push_back("Te");
        OxyCol.push_back("Po");
        
        std::vector<ID>::iterator iFind1 = std::find(OxyCol.begin(), OxyCol.end(), id1);
        std::vector<ID>::iterator iFind2 = std::find(OxyCol.begin(), OxyCol.end(), id2);
        
        if ((iFind1 != OxyCol.end() && bIdx1==2) && 
                (iFind2 != OxyCol.end() && bIdx2==3))
        {
            // two Oxy column atoms with sp3 orbiting
            //std::cout << "SetOneSP3OxyColumnBond(bIdx1, bIdx2, 2, 90.0)" << std::endl;
            //std::cout << "atom 1 "  << allAtoms[tIdx1].id 
            //          << " atom 2 " << allAtoms[tIdx2].id << std::endl;
            SetOneSP3OxyColumnBond(tIdx2, tIdx1, 3, 90.0);   
        }
        else if ((iFind1 != OxyCol.end() && bIdx1==2) && bIdx2==3)  
        {
            //std::cout << "SetOneSP3OxyColumnBond(bIdx1, bIdx2, 12, 180.0)" << std::endl;
            //std::cout << "atom 1 sp2 " << allAtoms[tIdx1].id 
            //        << " atom 2 sp3 " << allAtoms[tIdx2].id << std::endl;
            SetOneSP3OxyColumnBond(tIdx2, tIdx1, 3, 180.0, tF);
            
        }
        else if ((iFind2 != OxyCol.end() && bIdx2==2) && bIdx1==3)  
        {
            //std::cout << "SetOneSP3OxyColumnBond" << std::endl;
            //std::cout << "atom 1 sp3 " << allAtoms[tIdx1].id 
            //          << "atom 2 sp2 " << allAtoms[tIdx2].id << std::endl;
            SetOneSP3OxyColumnBond(tIdx1, tIdx2, 3, 180.0, tF);
        }
        else if (bIdx1==2 && bIdx2==2)
        {
            // std::cout << "SetOneSP2SP2Bond" << std::endl;
            SetOneSP2SP2Bond(tIdx1, tIdx2);
            
        }
        else if ((bIdx1==2 && bIdx2==3) || (bIdx1==3 && bIdx2==2))
        {
            if (bIdx1==2 && bIdx2==3)
            {
                //std::cout << "SetOneSP2SP3Bond(bIdx1, bIdx2)" << std::endl;
                //std::cout << " atom 1 sp2 " << allAtoms[tIdx1].id 
                //          << " atom 2 sp3 " << allAtoms[tIdx2].id << std::endl;
                SetOneSP2SP3Bond(tIdx1, tIdx2);
            }
            else
            {  
                //std::cout << "SetOneSP2SP3Bond(bIdx1, bIdx2)" << std::endl;
                //std::cout << " atom 1 sp3 " << allAtoms[tIdx1].id 
                //          << " atom 2 sp2 " << allAtoms[tIdx2].id << std::endl;
                SetOneSP2SP3Bond(tIdx2, tIdx1);
            }
        }
        else if (bIdx1==3 && bIdx2==3)
        {
            //std::cout << "SetOneSP3SP3Bond" << std::endl;
            
            SetOneSP3SP3Bond(tIdx1, tIdx2, tF);
        }  
    }
    
    void DictCifFile::setAllTorsions()
    {
        // loop over all bonds to get the torsion angles
        for (std::vector<BondDict>::iterator iABo= allBonds.begin();
                iABo != allBonds.end(); iABo++)
        {
            if ((int)iABo->fullAtoms.size() ==2)
            {
                std::vector<int> tPos;
                for(std::map<ID, int>::iterator iAM=iABo->fullAtoms.begin();
                        iAM !=iABo->fullAtoms.end(); iAM++)
                {
                    //std::cout << "atom " << iAM->first << std::endl;
                    // std::cout << " connected to " << (int)allAtoms[iAM->second].connAtoms.size()
                    //        << " atoms " << std::endl;
                    if ((int)allAtoms[iAM->second].connAtoms.size() > 1)
                    {
                        tPos.push_back(iAM->second);
                    }
                }
                // std::cout << "tPos.size() " << (int)tPos.size() << std::endl;
                if((int)tPos.size() ==2)
                {
                    setTorsionFromOneBond(tPos[0], tPos[1]);
                }
            }
        }
        
        std::cout << "All torsions have been setup " << std::endl;
        
    }
        
    void DictCifFile::setAllTorsions2()
    {
        
        std::vector<int>  tDone;
        
        // First set all torsion within all rings
        for (std::vector<RingDict>::iterator iR=allRingsV.begin();
                iR != allRingsV.end(); iR++)
        {
            setAllTorsionsInOneRing(tDone, *iR);
        }   
        
        // find all torsion not involved rings
        for (std::vector<BondDict>::iterator iABo= allBonds.begin();
                iABo != allBonds.end(); iABo++)
        {
            if (std::find(tDone.begin(), tDone.end(), iABo->seriNum)==tDone.end())
            {
                std::vector<int> tPos;
                for(std::map<ID, int>::iterator iAM=iABo->fullAtoms.begin();
                        iAM !=iABo->fullAtoms.end(); iAM++)
                {
                    //std::cout << "atom " << iAM->first << std::endl;
                    // std::cout << " connected to " << (int)allAtoms[iAM->second].connAtoms.size()
                    //        << " atoms " << std::endl;
                    if ((int)allAtoms[iAM->second].connAtoms.size() > 1)
                    {
                        tPos.push_back(iAM->second);
                    }
                }
                // std::cout << "tPos.size() " << (int)tPos.size() << std::endl;
                if((int)tPos.size() ==2)
                {
                    setTorsionFromOneBond(tPos[0], tPos[1]);
                }
            }
        }
        
        std::cout << "All torsions have been setup " << std::endl;
        
    }
    
    void DictCifFile::setAllTorsionsInOneRing(std::vector<int> & tBs, 
                                              RingDict         & tR)
    {
        
        std::vector<int> tAs, tLinkA, tBos;
        
        // A list of idx of atoms in the ring
        // std::cout << " atoms in the ring are: " << std::endl;
        for (int i=0; i < (int)tR.atoms.size(); i++)
        {
            tAs.push_back(tR.atoms[i].seriNum);
            // std::cout << tR.atoms[i].seriNum << std::endl;
        }
        
        tLinkA.push_back(tR.atoms[0].seriNum);
        
        int iCur   =tR.atoms[0].seriNum;
        int iLoop  =1;
        // std::cout << "ring rep " << tR.rep << " size " << (int)tAs.size() << std::endl;
       
        
        while ((int)tLinkA.size() < (int)tAs.size()
               && iLoop < (int)tAs.size())
        {   
            // std::cout << "atom " << allAtoms[iCur].seriNum << std::endl;
            for (std::vector<int>::iterator iC=allAtoms[iCur].connAtoms.begin();  
                    iC!=allAtoms[iCur].connAtoms.end(); iC++)
            {
                // std::cout << "connection " << *iC << std::endl;
                if(std::find(tAs.begin(), tAs.end(), *iC) !=tAs.end()
                   && std::find(tLinkA.begin(), tLinkA.end(), *iC) ==tLinkA.end())
                {
                    
                    int iB=getBond(allBonds, allAtoms[iCur].seriNum, allAtoms[*iC].seriNum);
                    if (iB >=0)
                    {
                        // std::cout << "find " << *iC << std::endl;
                        tBos.push_back(iB);
                        tLinkA.push_back(*iC);
                        iCur=*iC;
                        break;
                    }
                }
            }
            
            //double protection 
            iLoop++;
        }
        // Last bond
        /*
        int idxL = (int)tLinkA.size()-1;
        int iB=getBond(allBonds, allAtoms[tLinkA[0]].seriNum, 
                       allAtoms[tLinkA[idxL]].seriNum);
        if (iB >=0 && std::find(tBos.begin(), tBos.end(), iB)==tBos.end())
        {
            tBos.push_back(iB);
        }
        */
        if ((int)tLinkA.size() != (int)tR.atoms.size())
        {
            std::cout << "could not all linked atoms in ring " << tR.rep << std::endl;
            std::cout << "atoms in the ring are : " << std::endl;
            for (std::vector<AtomDict>::iterator iA=tR.atoms.begin(); 
                    iA !=tR.atoms.end(); iA++)
            {
                std::cout << iA->id << ",\t";
            }
            std::cout << std::endl;
            std::cout << "the links (bonds) found: " << std::endl;
            for (std::vector<int>::iterator iL=tLinkA.begin(); 
                    iL != tLinkA.end(); iL++)
            {
                std::cout << allAtoms[*iL].id << std::endl;
            }
        }
        //std::cout << "ring reps : " << tR.rep << std::endl;
        //std::cout << "it has " << (int)tBos.size() << std::endl;
        // Now we have bonds in the ring in sequence
        for (int i=0; i < (int)tBos.size(); i++)
        {
            std::string flip;
            if (i%2==0)
            {
                flip = "even";
            }
            else
            {
                flip = "odd";
            }
            
            std::vector<int> tPos;
            for(std::map<ID, int>::iterator iAM=allBonds[tBos[i]].fullAtoms.begin();
                        iAM !=allBonds[tBos[i]].fullAtoms.end(); iAM++)
            {
                //std::cout << "atom " << iAM->first << std::endl;
                // std::cout << " connected to " << (int)allAtoms[iAM->second].connAtoms.size()
                //        << " atoms " << std::endl;
                if ((int)allAtoms[iAM->second].connAtoms.size() > 1)
                {
                    tPos.push_back(iAM->second);
                }
            }
            // std::cout << "tPos.size() " << (int)tPos.size() << std::endl;
            if((int)tPos.size() ==2)
            {
                std::cout << "set torsion angles for the bond of atoms "
                        << allAtoms[tPos[0]].id << " and " 
                        << allAtoms[tPos[1]].id 
                        << "with " << flip << std::endl;
                        
                setTorsionFromOneBond(tPos[0], tPos[1], flip);
                tBs.push_back(tBos[i]);
            }
        }
    }
    
    void DictCifFile::setTorsionIdxFromOneBond(int tIdx1, int tIdx2)
    {
                
        std::cout << "For the bond consisting of atoms  " << allAtoms[tIdx1].id 
                  << " and " << allAtoms[tIdx2].id << std::endl
                  << "It has following torsion angles: " << std::endl; 
                
                
       for (std::vector<int>::iterator iAt1= allAtoms[tIdx1].connAtoms.begin();
            iAt1 != allAtoms[tIdx1].connAtoms.end(); iAt1++)
       {
           for (std::vector<int>::iterator iAt2 = allAtoms[tIdx2].connAtoms.begin();
                iAt2 != allAtoms[tIdx2].connAtoms.end(); iAt2++)
           {
               if (*iAt1 != tIdx2 && *iAt2 != tIdx1)
               {
                   TorsionDict aTorsion;         
                   aTorsion.atoms.push_back(*iAt1);
                   aTorsion.atoms.push_back(tIdx1);
                   aTorsion.atoms.push_back(tIdx2);
                   aTorsion.atoms.push_back(*iAt2);
                   std::cout << "Torsion: " << allAtoms[*iAt1].id 
                             << ", " << allAtoms[tIdx1].id 
                             << ", " << allAtoms[tIdx2].id
                             << ", " << allAtoms[*iAt2].id << std::endl; 
                                    
                   allTorsions.push_back(aTorsion);
               }
           }
        }
        
    }
    
    void DictCifFile::getChiralInfo(std::vector<std::string> tF)
    {
        // check which entries exist in an chiral section 
        
        
        if ((int)tF.size() ==1 &&
                tF[0].find("_chem_comp_chir.")!=std::string::npos)
        {
            std::vector<std::string> tF1;
            tF1 = StrTokenize(TrimSpaces(tF[0]), '.');
            if((int)tF1.size() ==2)
            {
               if(tF1[1].find("comp_id") !=std::string::npos)
               {
                   //std::cout << "resName" << std::endl;
                   hasProps["chiral"].insert(std::pair<std::string, int>("resName",curBlockLine));
                   //std::cout << (int) hasProps["chiral"].size() << std::endl;
                   curBlockLine++;
               } 
               else if(tF1[1].compare("id")==0)
               {
                   //std::cout << "id" << std::endl;
                   hasProps["chiral"].insert(std::pair<std::string, int>("id",curBlockLine));
                   //std::cout << (int) hasProps["chiral"].size() << std::endl;
                   curBlockLine++;
               }
               else if(tF1[1].find("atom_id_centre") !=std::string::npos)
               {
                   hasProps["chiral"].insert(std::pair<std::string, int>("archID",curBlockLine));
                   //std::cout << (int) hasProps["chiral"].size() << std::endl;
                   curBlockLine++;
               }
               else if(tF1[1].find("atom_id_1") !=std::string::npos)
               {
                   hasProps["chiral"].insert(std::pair<std::string, int>("atom1",curBlockLine));
                   //std::cout << (int) hasProps["chiral"].size() << std::endl;
                   curBlockLine++;
               }
               else if(tF1[1].find("atom_id_2") !=std::string::npos)
               {
                   hasProps["chiral"].insert(std::pair<std::string, int>("atom2",curBlockLine));
                   //std::cout << (int) hasProps["chiral"].size() << std::endl;
                   curBlockLine++;
               }
               else if(tF1[1].find("atom_id_3") !=std::string::npos)
               {
                   hasProps["chiral"].insert(std::pair<std::string, int>("atom3",curBlockLine));
                   //std::cout << (int) hasProps["chiral"].size() << std::endl;
                   curBlockLine++;
               }
               else if(tF1[1].find("volume_sign") !=std::string::npos)
               {                    
                   
                   hasProps["chiral"].insert(std::pair<std::string, int>("volume_sign",curBlockLine));
                   //std::cout << (int) hasProps["chiral"].size() << std::endl;
                   curBlockLine++;
               }
            }
        }
        
            // Now the content lines with entries match to hasProps["chiral"]
            // std::cout << (int)tF.size() << std::endl;
            if ((int)tF.size() == (int)hasProps["chiral"].size() 
                && (int)tF.size() >2 && tF[0].find("#") ==std::string::npos)
            {
                itsCurChiral = new ChiralDict();
                  
                if (hasProps["chiral"].find("id") != hasProps["chiral"].end())
                {
                    cleanSymbol(tF[hasProps["chiral"]["id"]], "\"");
                    itsCurChiral->id = TrimSpaces(tF[hasProps["chiral"]["id"]]);
                }
                if (hasProps["chiral"].find("archID") != hasProps["chiral"].end())
                {
                    cleanSymbol(tF[hasProps["chiral"]["archID"]], "\"");
                    ID tAtId = TrimSpaces(tF[hasProps["chiral"]["archID"]]);
                    itsCurChiral->archID = tAtId;
                    int iPos = atomPosition(tAtId);
                    if (iPos !=-1)
                    {
                        itsCurChiral->archPos = iPos;
                        itsCurChiral->atoms.push_back(iPos);
                        allAtoms[iPos].inChirals.push_back(itsCurChiralSeriNum);
                    }
                    else
                    {
                        std::cout << "Chiral definition error? Could not find " 
                                  << tAtId  << " existing atoms " << std::endl;
                        exit(1);
                    }    
                }
                if (hasProps["chiral"].find("atom1") != hasProps["chiral"].end())
                {
                    cleanSymbol(tF[hasProps["chiral"]["atom1"]], "\"");
                    ID tAtId = TrimSpaces(tF[hasProps["chiral"]["atom1"]]);
                    int iPos = atomPosition(tAtId);
                    if(iPos !=-1)
                    {
                        itsCurChiral->atoms.push_back(iPos);  
                    }
                    else
                    {
                        std::cout << "Chiral definition error? Could not find " 
                                   << tAtId  << " existing atoms " << std::endl;
                         exit(1);
                    }
                }
                if (hasProps["chiral"].find("atom2") != hasProps["chiral"].end())
                {
                    cleanSymbol(tF[hasProps["chiral"]["atom2"]], "\"");
                    ID tAtId = TrimSpaces(tF[hasProps["chiral"]["atom2"]]);
                    int iPos = atomPosition(tAtId);
                    if(iPos !=-1)
                    {
                        itsCurChiral->atoms.push_back(iPos);  
                    }
                    else
                    {
                        std::cout << "Chiral definition error? Could not find " 
                                  << tAtId  << " existing atoms " << std::endl;
                           exit(1);
                    }
                }
                    
                if (hasProps["chiral"].find("atom3") != hasProps["chiral"].end())
                {
                    cleanSymbol(tF[hasProps["chiral"]["atom3"]], "\"");
                    ID tAtId = TrimSpaces(tF[hasProps["chiral"]["atom3"]]);
                    int iPos = atomPosition(tAtId);
                    if(iPos !=-1)
                    {
                        itsCurChiral->atoms.push_back(iPos);  
                    }
                    else
                    {
                        std::cout << "Chiral definition error? Could not find " 
                                  << tAtId  << " existing atoms " << std::endl;
                        exit(1);
                    }
                }
                if (hasProps["chiral"].find("volume_sign") != hasProps["chiral"].end())
                {
                    cleanSymbol(tF[hasProps["chiral"]["volume_sign"]], "\"");
                    itsCurChiral->sign = TrimSpaces(tF[hasProps["chiral"]["volume_sign"]]);
                    if (itsCurChiral->sign =="positiv")
                    {
                        allAtoms[itsCurChiral->atoms[0]].chiralIdx = 1;
                    }
                    else if (itsCurChiral->sign =="negativ")
                    {
                        allAtoms[itsCurChiral->atoms[0]].chiralIdx = -1;
                    }
                    else if (itsCurChiral->sign =="both")
                    {
                        allAtoms[itsCurChiral->atoms[0]].chiralIdx = 2;
                    }
                    else
                    {
                        allAtoms[itsCurChiral->atoms[0]].chiralIdx = 0;
                    }  
                }
                
                   
                if((int)itsCurChiral->atoms.size() ==4)
                {  
                    // Check which atom is not listed in the chiral list
                    int idxA=-1;
                    std::vector<int>::iterator iFA;
                    for (std::vector<int>::iterator iTA=allAtoms[itsCurChiral->atoms[0]].connAtoms.begin();
                            iTA != allAtoms[itsCurChiral->atoms[0]].connAtoms.end(); iTA++)
                    {
                        iFA=std::find(itsCurChiral->atoms.begin(), itsCurChiral->atoms.end(), *iTA);
                        if (iFA ==itsCurChiral->atoms.end())
                        {
                            idxA=*iTA;
                            break;
                        }
                    }
                    
                    if (allAtoms[itsCurChiral->atoms[0]].chiralIdx==1 ||
                        allAtoms[itsCurChiral->atoms[0]].chiralIdx==2)
                    {
                        itsCurChiral->mutTable[itsCurChiral->atoms[1]].push_back(itsCurChiral->atoms[3]);
                        itsCurChiral->mutTable[itsCurChiral->atoms[1]].push_back(itsCurChiral->atoms[2]);
                       
                        itsCurChiral->mutTable[itsCurChiral->atoms[2]].push_back(itsCurChiral->atoms[1]);
                        itsCurChiral->mutTable[itsCurChiral->atoms[2]].push_back(itsCurChiral->atoms[3]);
                        
                        itsCurChiral->mutTable[itsCurChiral->atoms[3]].push_back(itsCurChiral->atoms[2]);
                        itsCurChiral->mutTable[itsCurChiral->atoms[3]].push_back(itsCurChiral->atoms[1]);
                        if (idxA !=-1)
                        {
                            itsCurChiral->mutTable[itsCurChiral->atoms[1]].push_back(idxA);
                            itsCurChiral->mutTable[itsCurChiral->atoms[2]].push_back(idxA);
                            itsCurChiral->mutTable[itsCurChiral->atoms[3]].push_back(idxA);
                            
                            // table for the atom not in the chiral list
                            for(int i=1; i < (int)itsCurChiral->atoms.size(); i++)
                            {
                                itsCurChiral->mutTable[idxA].push_back(itsCurChiral->atoms[i]);
                            }
                        }         
                    }
                    
                    if (allAtoms[itsCurChiral->atoms[0]].chiralIdx==-1)
                    {
                        itsCurChiral->mutTable[itsCurChiral->atoms[1]].push_back(itsCurChiral->atoms[2]);
                        itsCurChiral->mutTable[itsCurChiral->atoms[1]].push_back(itsCurChiral->atoms[3]);
                       
                        itsCurChiral->mutTable[itsCurChiral->atoms[2]].push_back(itsCurChiral->atoms[3]);
                        itsCurChiral->mutTable[itsCurChiral->atoms[2]].push_back(itsCurChiral->atoms[1]);
                        
                        itsCurChiral->mutTable[itsCurChiral->atoms[3]].push_back(itsCurChiral->atoms[1]);
                        itsCurChiral->mutTable[itsCurChiral->atoms[3]].push_back(itsCurChiral->atoms[2]);
                        if (idxA !=-1)
                        {
                            itsCurChiral->mutTable[itsCurChiral->atoms[1]].push_back(idxA);
                            itsCurChiral->mutTable[itsCurChiral->atoms[2]].push_back(idxA);
                            itsCurChiral->mutTable[itsCurChiral->atoms[3]].push_back(idxA);
                            // table for atom not in the chiral list
                           for(int i=(int)itsCurChiral->atoms.size()-1; i>0; i--)
                           {
                               itsCurChiral->mutTable[idxA].push_back(itsCurChiral->atoms[i]);
                           }
                        }
                    }
                    
                    // std::cout << "atom " <<  allAtoms[itsCurChiral->atoms[0]].id << " chiral sign is " 
                    //          << allAtoms[itsCurChiral->atoms[0]].chiralIdx << std::endl
                    //          << " mutable size " << (int)itsCurChiral->mutTable.size()
                    //          << std::endl;
                    
                    /*
                    for (std::map<int, std::vector<int> >::iterator iM=itsCurChiral->mutTable.begin();
                            iM !=itsCurChiral->mutTable.end(); iM++)
                    {
                        std::cout << allAtoms[iM->first].id << " mutable are: " << std::endl;
                        for (std::vector<int>::iterator iA=iM->second.begin(); 
                                iA != iM->second.end(); iA++)
                        {
                            std::cout << "atom ID " << allAtoms[*iA].id << "  serial number " << allAtoms[*iA].seriNum
                                      << std::endl;
                        }
                    }
                     */
                }
                
                allChirals.push_back(*itsCurChiral);
                itsCurChiralSeriNum++;
                //std::cout << "Number of chirals in the system " << (int)allChirals.size() << std::endl;
                delete itsCurChiral;
                itsCurChiral = NULL;
            }
        
    }
    
    void DictCifFile::setDefaultCoordGeos()
    {
        DefaultCoordGeos[2]  = "LINEAR";
        DefaultCoordGeos[3]  = "TRIGONAL-PLANAR";
        DefaultCoordGeos[4]  = "TETRAHEDRAL";
        DefaultCoordGeos[5]  = "TRIGONAL-BIPYRAMID";
        DefaultCoordGeos[6]  = "OCTAHEDRAL";
        DefaultCoordGeos[7]  = "CAPPED-OCTAHEDRAL";
        DefaultCoordGeos[8]  = "CUBIC";
        DefaultCoordGeos[9]  = "TRICAPPED-TRIGONAL-PRISMATIC";
        DefaultCoordGeos[10] = "BICAPPED-SQUARE-ANTIPRISMATIC";
        DefaultCoordGeos[11] = "ALL-FACE-CAPPED-TRIGONAL-PRISMATIC";
        DefaultCoordGeos[12] = "CUBOCTAHEDRON";
        
        //std::string clibMonDir(std::getenv("CLIBD_MON"));
        //std::string metDefCoordGeoFileName = clibMonDir + "/allMetalDefCoordGeos.table";
        std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
        std::string metDefCoordGeoFileName = clibMonDir + "/tables/allMetalDefCoordGeos.table";
        std::ifstream metDefCoordGeoFile(metDefCoordGeoFileName.c_str());
        
        if(metDefCoordGeoFile.is_open())
        {
            std::string tRecord="";
            
            while(!metDefCoordGeoFile.eof())
            {
                std::getline(metDefCoordGeoFile, tRecord);
                tRecord = TrimSpaces(tRecord);
                std::vector<std::string> tBuf;
                StrTokenize(tRecord, tBuf);
                
                if ((int)tBuf.size() ==3)
                {
                    int cn = StrToInt(TrimSpaces(tBuf[1]));
                    
                    DefaultCoordGeos2[TrimSpaces(tBuf[0])][cn]=TrimSpaces(tBuf[2]);
                }
            }
            metDefCoordGeoFile.close();
        }
        
    }   
    
    void DictCifFile::getMetalCNGeo(std::vector<std::string> tF)
    {
        // check which entries exist in an metal atom
        if ((int)tF.size() ==1 && 
            tF[0].find("_chem_comp_metal.")!=std::string::npos)
        {
            std::vector<std::string> tF1;
            tF1 = StrTokenize(TrimSpaces(tF[0]), '.');
            if((int)tF1.size() ==2)
            {
                if(tF1[1].find("comp_id") !=std::string::npos)
                {
                    hasProps["metal"].insert(std::pair<std::string, int>("resName",curBlockLine));
                    //std::cout << curBlockLine << std::endl;
                    curBlockLine++;
                }
                else if(tF1[1].find("atom_id") !=std::string::npos)
                {
                    hasProps["metal"].insert(std::pair<std::string, int>("id",curBlockLine));
                    //std::cout << curBlockLine << std::endl;
                    curBlockLine++;
                }
                else if(tF1[1].find("coord_num") !=std::string::npos)
                {
                    hasProps["metal"].insert(std::pair<std::string, int>("coordNum",curBlockLine));
                    //std::cout << curBlockLine << std::endl;
                    curBlockLine++;
                }
                else if(tF1[1].find("coord_geo") !=std::string::npos)
                {
                    hasProps["metal"].insert(std::pair<std::string, int>("coordGeo",curBlockLine));
                    //std::cout << curBlockLine << std::endl;
                    curBlockLine++;
                }
            }
        }
        if ((int)tF.size() == (int)hasProps["metal"].size()
            && (int)tF.size() > 2 && tF[0].find("#") ==std::string::npos)
        {
            if (hasProps["metal"].find("id") != hasProps["metal"].end())
            {
                cleanSymbol(tF[hasProps["metal"]["id"]], "\"");
                int tAtomIdx = atomPosition(tF[hasProps["metal"]["id"]]);
                if (tAtomIdx !=-1)
                {
                    allAtoms[tAtomIdx].isMetal = true;
                    cleanSymbol(tF[hasProps["metal"]["coordGeo"]], "\"");
                    allAtoms[tAtomIdx].metalGeo = tF[hasProps["metal"]["coordGeo"]];
                    cleanSymbol(tF[hasProps["metal"]["coordNum"]], "\"");
                    int iCN = StrToInt(tF[hasProps["metal"]["coordNum"]]);
                    if (iCN !=(int)allAtoms[tAtomIdx].connAtoms.size())
                    {
                        std::cout << "metal Atom " << allAtoms[tAtomIdx].id 
                                  << " is bonding to " << (int)allAtoms[tAtomIdx].connAtoms.size()
                                  << " atoms. " <<  std::endl 
                                  << "But the coordination number defined in the input is "
                                  << iCN << " Please check the input file. The job process stops"
                                  << std::endl;
                        exit(1);
                    }
                    std::cout << "Atom " << allAtoms[tAtomIdx].id 
                              << " is metal. its coordination number is " 
                              << iCN  << std::endl
                              << " and its input coordination geometry is "
                              << allAtoms[tAtomIdx].metalGeo << std::endl;
                }
            }
        }
    }
        
    
    short DictCifFile::transOrder(std::string tO)
    {
        std::string tS = TrimSpaces(tO);
        short tOrder = 0;
        StrLower(tS);
        if (tS.compare("single") == 0 || tS.compare(0,4,"sing") == 0 )
        {
            tOrder = 1;
        }
        else if (tS.compare("double") == 0 || tS.compare(0,4,"doub") == 0 )
        {
            tOrder = 2;
        }
        else if (tS.compare("triple") == 0 || tS.compare(0,4,"trip") == 0 )
        {
            tOrder = 3;
        }
        
        // and tOrder could be 4 
        // deloc = 0 as well
        
        return tOrder;
        
    }
    
    
    void DictCifFile::outAtomBloc()
    {
        
    }
    
    void DictCifFile::outSystem()
    {
        
    }
    
    void DictCifFile::transCoordsPdbToCif(DictPDBFile       & tPdbObj)
    {
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA !=allAtoms.end(); iA++)
        {
            if (tPdbObj.allHetAtmList.find(iA->id) !=tPdbObj.allHetAtmList.end())
            {
                if ((tPdbObj.allHetAtmList[iA->id].size() ==iA->coords.size())
                     && iA->coords.size() !=0)
                {
                    for (unsigned i=0; i < iA->coords.size(); iA++)
                    {
                        iA->coords[i] = tPdbObj.allHetAtmList[iA->id][i];
                    }
                }
            }
            
        }
        
    }
    
    extern void outMMCif(FileName tFName, 
                         ID tMonoRootName,
                         ChemComp  &         tPropComp,
                         std::vector<AtomDict>& tAtoms,
                         // std::vector<int>    & tHydroAtoms,
                         std::vector<BondDict>& tBonds, 
                         std::vector<AngleDict>& tAngs, 
                         std::vector<TorsionDict>& tTorsions, 
                         std::vector<RingDict> & tRings, 
                         std::vector<PlaneDict>& tPlas, 
                         std::vector<ChiralDict>& tChs,
                         const   double           tUBS,
                         const   double           tLBS,
                         const   double           tUAS,
                         const   double           tLAS,
                         std::map<int, std::map<std::string,
                         std::map<std::string, double > > > 
                         &  tHDistMap)
    {
        
        /*
        for (std::vector<AtomDict>::iterator iAt= tAtoms.begin();
                iAt!=tAtoms.end(); iAt++)
        {
            std::cout << "atom " << iAt->id << " has charge of "
                      << iAt->formalCharge << std::endl;
        }
        
        
        // newly added 
        if(tAtoms.size()> 0 && tBonds.size() > 0)
        {
            //std::cout << "Kekulize the system " << std::endl;
            KekulizeMol aKTool;
            aKTool.execute(tAtoms, 
                           tBonds,
                           tRings);
            //std::cout << "Kekulize done " << std::endl;
        }
        */

        for (std::vector<AtomDict>::iterator iA=tAtoms.begin();
                iA !=tAtoms.end(); iA++)
        {
            // std::cout << iA->id << std::endl;
            if (iA->id.find("\'") !=std::string::npos)
            {
                iA->id = "\"" + iA->id + "\"";
            }
            // std::cout << iA->id << std::endl;
        }
        /*
        for (std::vector<AtomDict>::iterator iAt= tAtoms.begin();
                iAt!=tAtoms.end(); iAt++)
        {
            std::cout << "atom " << iAt->id << " has charge of "
                      << iAt->formalCharge << std::endl;
        }
        */
        
        std::vector<std::string> aSetStrs;
        StrTokenize(tFName, aSetStrs, '.');
        
        
        std::ofstream outRestrF(tFName);
        
        if(outRestrF.is_open())
        {
            
            srand((unsigned)std::time( NULL ));
            // Temp 
            // 1. Global section 
            outRestrF << "global_" << std::endl
                    << "_lib_name         ?" << std::endl
                    << "_lib_version      ?" << std::endl
                    << "_lib_update       ?" << std::endl;
        
            
            
            // 'LIST OF MONOMERS' section
            
            
            std::string longName =tMonoRootName.substr(0,3);
            std::string sName =tMonoRootName.substr(0,3);
            
            //StrUpper(longName);
            
            
            ID ligType = "non-polymer";
            
            for (std::vector<LIBMOL::RingDict>::iterator iR=tRings.begin();
                    iR != tRings.end(); iR++)
            {   
                if (iR->sugarType.compare("pyranose")==0)
                {
                    ligType = "pyranose";
                    break;
                }
             
            }
            
            std::vector<ID>  aAATab;
            initAminoAcidTab(aAATab);
            if (isAminoAcid(aAATab, longName) && longName.find("PRO")==std::string::npos)
            {
                ligType = "L-peptide";
            }
            
   
            int nH = getHAtomNum(tAtoms);
            
            
            outRestrF << "# ------------------------------------------------" << std::endl
                    << "#" << std::endl
                    << "# ---   LIST OF MONOMERS ---" << std::endl
                    << "#" << std::endl
                    << "data_comp_list" << std::endl
                    << "loop_" << std::endl
                    << "_chem_comp.id" << std::endl
                    << "_chem_comp.three_letter_code" << std::endl
                    << "_chem_comp.name" << std::endl
                    << "_chem_comp.group" << std::endl
                    << "_chem_comp.number_atoms_all" << std::endl
                    << "_chem_comp.number_atoms_nh"  << std::endl
                    << "_chem_comp.desc_level" << std::endl;
            
            //if (tPropComp.id !=NullString)
            //{
            //    outRestrF << tPropComp.id  <<"\t"<< tPropComp.code << "\t" 
            //              << tPropComp.name << "\t" << tPropComp.group << "\t" 
            //              << tPropComp.numAtoms << "\t" 
            //              << tPropComp.numH << "\t."
                          // << (int)tAtoms.size()-(int)tHydroAtoms.size() << "\t."
            //              << std::endl;
                
            //}
            //else
            //{
            
            
            outRestrF << longName <<"\t"<< sName << "\t" << "'.\t\t'\t"
                      << ligType << "\t" << (int)tAtoms.size() << "\t" 
                      << (int)tAtoms.size()- nH << "\t."
                          // << (int)tAtoms.size()-(int)tHydroAtoms.size() << "\t."
                      << std::endl;
            
            outRestrF <<"# ------------------------------------------------------" << std::endl
                      <<"# ------------------------------------------------------" << std::endl
                      <<"#" << std::endl
                      <<"# --- DESCRIPTION OF MONOMERS ---" << std::endl
                      <<"#" << std::endl
                      <<"data_comp_" << longName << std::endl
                      <<"#" << std::endl; 
        
            if (tAtoms.size() >0)
            {
                // atom info section           
                outRestrF << "loop_" << std::endl
                          << "_chem_comp_atom.comp_id" << std::endl
                          << "_chem_comp_atom.atom_id" << std::endl
                          << "_chem_comp_atom.type_symbol" << std::endl
                          << "_chem_comp_atom.type_energy" << std::endl
                          << "_chem_comp_atom.charge" << std::endl
                          << "_chem_comp_atom.x" << std::endl
                          << "_chem_comp_atom.y" << std::endl
                          << "_chem_comp_atom.z" << std::endl;
                
                
                
                for (std::vector<AtomDict>::iterator iA = tAtoms.begin();
                        iA != tAtoms.end(); iA++)
                {
                    //double r1 =  (double) rand()/RAND_MAX;
                    //double r2 =  (double) rand()/RAND_MAX;
                    //double r3 =  (double) rand()/RAND_MAX;
                    REAL tCharge =0.0;
                    //if (iA->charge !=0.0)
                    //{
                    //    tCharge = iA->charge;
                    //}
                    if (iA->formalCharge !=0.0)
                    {
                        tCharge = iA->formalCharge;
                    }
                    
                    std::string strCharge = TrimSpaces(RealToStr(tCharge));
                    if (strCharge.find(".") !=strCharge.npos)
                    {
                        std::vector<std::string> tVec;
                        StrTokenize(strCharge, tVec, '.');
                        strCharge = tVec[0];
                    }
                    
                    StrUpper(iA->chemType);
                    StrUpper(iA->ccp4Type);
                    
                    outRestrF << longName
                              << std::setw(12) << iA->id 
                              << std::setw(6) << iA->chemType 
                              << std::setw(6) << iA->ccp4Type 
                              << std::setw(8) << strCharge 
                              << std::setw(12) << std::setprecision(3) << std::fixed 
                              << iA->coords[0] 
                              << std::setw(12) << std::setprecision(3) << std::fixed 
                              << iA->coords[1]
                              << std::setw(12) << std::setprecision(3) << std::fixed 
                              << iA->coords[2] << std::endl;
                    
                }
            }
            
            if (tBonds.size() >0)
            {   
                // Bond sections 
                outRestrF << "loop_" << std::endl
                          << "_chem_comp_bond.comp_id" << std::endl
                          << "_chem_comp_bond.atom_id_1" << std::endl
                          << "_chem_comp_bond.atom_id_2" << std::endl
                          << "_chem_comp_bond.type" << std::endl
                          << "_chem_comp_bond.aromatic"  << std::endl
                          << "_chem_comp_bond.value_dist_nucleus" << std::endl
                          << "_chem_comp_bond.value_dist_nucleus_esd" << std::endl
                          << "_chem_comp_bond.value_dist"<< std::endl
                          << "_chem_comp_bond.value_dist_esd" << std::endl;
                          // << "_chem_comp_bond.exact_cod_dist" << std::endl;
               
                
                for (std::vector<BondDict>::iterator iB=tBonds.begin();
                          iB !=tBonds.end(); iB++)
                {
                    std::string tAr;
                    
                    
                    StrLower(iB->order);
                    
                    if (iB->order.find("arom") !=iB->order.npos)
                    {
                        tAr       = "y";
                    }
                    else
                    {
                        unifyStrForOrder(iB->order);
                        tAr = "n";
                    }
                    
                    std::string outOrder;
                    if (iB->orderNI.size() !=0)
                    {
                        outOrder = iB->orderNI;
                    }
                    else
                    {
                        outOrder = iB->order;
                    }
                    
                    double aSigV;
                    if (iB->sigValue < tLBS)
                    {
                        aSigV = tLBS;
                    }
                    else if (iB->sigValue > tUBS)
                    {
                        aSigV = tUBS;
                    }
                    else
                    {
                        aSigV  = iB->sigValue;
                    }
                    
                    // Add Proton-X distances if they exist
                    int idxH = -1, idxX=-1;
                    
                    if (tAtoms[iB->atomsIdx[0]].chemType.compare("H")==0)
                    {
                        idxH = iB->atomsIdx[0];
                        idxX = iB->atomsIdx[1];
                    }
                    else if (tAtoms[iB->atomsIdx[1]].chemType.compare("H")==0)
                    {
                        idxH = iB->atomsIdx[1];
                        idxX = iB->atomsIdx[0];
                    }
                    
                    double aProtD      = 0.0;
                    double aProtD_siga = 0.0;
                            
                               
                    if (idxH > -1)
                    {
                        
                        if (tHDistMap[2].find(tAtoms[idxH].formType[1])
                             !=tHDistMap[2].end())
                        {
                            aProtD = tHDistMap[2][tAtoms[idxH].formType[1]]["pDistNeu"];
                            aProtD_siga 
                            = tHDistMap[2][tAtoms[idxH].formType[1]]["pDistSigaNeu"];
                        }
                        else if (tHDistMap[1].find(tAtoms[idxX].chemType)
                             !=tHDistMap[1].end())
                        {
                            aProtD = tHDistMap[1][tAtoms[idxX].chemType]["pDist"];
                            aProtD_siga 
                            = tHDistMap[1][tAtoms[idxX].chemType]["pDistSiga"];
                        }
                        else
                        {
                            aProtD = iB->value;
                            aProtD_siga = aSigV;
                        }
                    }
                    else
                    {
                        aProtD = iB->value;
                        aProtD_siga = aSigV;
                    }
                    
                    outRestrF <<  longName
                              << std::setw(12)  << tAtoms[iB->atomsIdx[0]].id  
                              << std::setw(12)  << tAtoms[iB->atomsIdx[1]].id  
                              //<< std::setw(12)  << iB->orderNK 
                              << std::setw(12)  << outOrder
                              << std::setw(8)   << tAr
                              << std::setw(10)  << std::setprecision(3)
                              << aProtD
                              << std::setw(8) << std::setprecision(4)
                              << aProtD_siga
                              << std::setw(10)  << std::setprecision(3)
                              << iB->value 
                              << std::setw(8) << std::setprecision(4)
                              << aSigV << std::endl;
                    
                    
                }
            }
           
            if (tAngs.size() > 0)
            {
                // Angle section
                outRestrF << "loop_" << std::endl
                          << "_chem_comp_angle.comp_id"   << std::endl
                          << "_chem_comp_angle.atom_id_1" << std::endl
                          << "_chem_comp_angle.atom_id_2" << std::endl
                          << "_chem_comp_angle.atom_id_3" << std::endl
                          << "_chem_comp_angle.value_angle"     << std::endl
                          << "_chem_comp_angle.value_angle_esd" << std::endl;
                          //  << "_chem_comp_angle.exact_cod_dist"  << std::endl;
                
                for (std::vector<AngleDict>::iterator iA=tAngs.begin();
                          iA != tAngs.end(); iA++)
                {
                    //for (std::vector<int>::iterator iAt=iA->atoms.begin();
                    //        iAt !=iA->atoms.end(); iAt++)
                    //{
                    // difference in comp atom definitions between cod and
                    // dictionary: inner-out1-out2(cod),
                    // atom1-atom2(center)-atom3(dictionary)
                    double aSigAV;
                    
                    if (iA->sigValue < tLAS)
                    {
                        aSigAV = tLAS;
                    }
                    else if (iA->sigValue > tUAS)
                    {
                        aSigAV = tUAS;
                    }
                    else
                    {
                        aSigAV = iA->sigValue;
                    }
                    
                    if (tAtoms[iA->atoms[0]].isMetal)
                    {
                        for (std::vector<REAL>::iterator iCA=iA->codAngleValues.begin();
                                iCA !=iA->codAngleValues.end(); iCA++)
                        {
                            outRestrF << tMonoRootName.substr(0,3) 
                                      << std::setw(12)
                                      << tAtoms[iA->atoms[1]].id 
                                      << std::setw(12)
                                      << tAtoms[iA->atoms[0]].id 
                                      << std::setw(12)
                                      << tAtoms[iA->atoms[2]].id 
                                      << std::setw(12) 
                                      << std::setprecision(3) 
                                      <<  *iCA << "    "
                                      << std::setw(8) << std::setprecision(2) 
                                      << aSigAV << std::endl;
                        }
                    }
                    else
                    {
                        outRestrF << longName 
                                  << std::setw(12) << tAtoms[iA->atoms[1]].id
                                  << std::setw(12) << tAtoms[iA->atoms[0]].id 
                                  << std::setw(12) << tAtoms[iA->atoms[2]].id;
                        outRestrF << std::setw(12) << std::setprecision(3) <<  iA->value
                                  << std::setw(8) << std::setprecision(2) << aSigAV 
                                  << std::endl;
                    }
                    /*
                    if(iA->hasCodValue)
                     {
                       outRestrF << "Yes " << std::endl;
                     }
                     else
                     {
                       outRestrF << "No "  << std::endl;
                     }
                     */
                }
            }
        
            // Torsion section 
            if((int)tTorsions.size() !=0)
            {
                outRestrF << "loop_" << std::endl
                          << "_chem_comp_tor.comp_id"         << std::endl
                          << "_chem_comp_tor.id"              << std::endl
                          << "_chem_comp_tor.atom_id_1"       << std::endl
                          << "_chem_comp_tor.atom_id_2"       << std::endl
                          << "_chem_comp_tor.atom_id_3"       << std::endl
                          << "_chem_comp_tor.atom_id_4"       << std::endl
                          << "_chem_comp_tor.value_angle"     << std::endl
                          << "_chem_comp_tor.value_angle_esd" << std::endl
                          << "_chem_comp_tor.period"          << std::endl;
            
                int idxTor = 1;
                //std::cout << "number of torsions " << (int)allTorsions.size() << std::endl;
            
                for (std::vector<TorsionDict>::iterator iT=tTorsions.begin();
                        iT !=tTorsions.end(); iT++)
                {
                    //std::string idxTorStr=IntToStr(idxTor);
                    //idxTorStr = "tor_" + idxTorStr;
                    // std::cout << "Torsion angle " << idxTor 
                    //          << " It contains " << (int)iT->atoms.size() << std::endl;
                     /*
                    std::cout << "atom sp in a torsion " << std::endl;
                    
                    std::cout << "atom " << tAtoms[iT->atoms[1]].id
                              << " : " << tAtoms[iT->atoms[1]].hybrid
                              << std::endl
                              << "atom " << tAtoms[iT->atoms[2]].id
                              << " : " << tAtoms[iT->atoms[1]].hybrid
                              << std::endl;
                    */
                    
                    std::string aTorSiga;
                    if (iT->id.find("const") !=std::string::npos
                        || iT->id.find("CONST") !=std::string::npos)
                    {
                        aTorSiga = "0.0";
                    }
                    else if (iT->id.find("sp2_sp2") !=std::string::npos)
                    {
                        aTorSiga = "5.0";
                    }
                    else if (tAtoms[iT->atoms[1]].hybrid.find("SP2") 
                             !=std::string::npos
                            && tAtoms[iT->atoms[2]].hybrid.find("SP2")
                             !=std::string::npos)
                    {
                        aTorSiga = "5.0";
                    }
                    else
                    {
                        aTorSiga = "10.0";
                    }
                    outRestrF << longName 
                              << std::setw(22) << iT->id
                              << std::setw(12)  << tAtoms[iT->atoms[0]].id 
                              << std::setw(12)  << tAtoms[iT->atoms[1]].id 
                              << std::setw(12)  << tAtoms[iT->atoms[2]].id 
                              << std::setw(12)  << tAtoms[iT->atoms[3]].id 
                              << std::setw(12) << std::setprecision(3) << iT->value  
                              << std::setw(8)  << aTorSiga
                              << std::setw(6)  << iT->period << std::endl;
                    idxTor++;        
                }
                
            }
            
            //  For chiral centers
            //bool l_ch = false;
            //if ((int)tChs.size() !=0)
            //{
            //    l_ch = true;
            //}
            //else
            //{
            //    for (int i_ch =0; i_ch < (int)tAtoms.size(); i_ch++)
            //    {
            //        if (tAtoms[i_ch].chiralIdx  > 0)
            //        {
            //            l_ch = true;
            //            break;
            //        }
            //    }
            //}
            
            //if (l_ch)
            if(tChs.size() !=0)
            {
                outRestrF << "loop_" << std::endl
                          << "_chem_comp_chir.comp_id" << std::endl
                          << "_chem_comp_chir.id" << std::endl
                          << "_chem_comp_chir.atom_id_centre" << std::endl
                          << "_chem_comp_chir.atom_id_1" << std::endl
                          << "_chem_comp_chir.atom_id_2" << std::endl
                          << "_chem_comp_chir.atom_id_3" << std::endl
                          << "_chem_comp_chir.volume_sign" << std::endl;
                
                // First the input chirals
                std::vector<ID>   inputChiralID;
                for (std::vector<ChiralDict>::iterator iCh = tChs.begin();
                        iCh != tChs.end(); iCh++)
                {
                   
                        inputChiralID.push_back(iCh->archID);
                        outRestrF << longName << "    " 
                                  << iCh->id  << "    ";
                        int numCh=0;
                        for (std::vector<int>::iterator iAt=iCh->atoms.begin();
                               iAt != iCh->atoms.end(); iAt++)
                        {
                            if (numCh < 4)
                            {
                                outRestrF << tAtoms[*iAt].id << "    ";
                                numCh++;
                            }
                        }
                        outRestrF << iCh->sign << std::endl;
                   
                }
                // New chiral that are not in the input list 
                /*
                int idxC =(int)inputChiralID.size();
                for (std::vector<AtomDict>::iterator iA=tAtoms.begin();
                          iA !=tAtoms.end(); iA++)
                {
                 
                    if (iA->chiralIdx == 1)
                    {
                        std::vector<ID>::iterator tFind;
                        tFind = std::find(inputChiralID.begin(), inputChiralID.end(), iA->id); 
                        if (tFind ==inputChiralID.end() && iA->chemType !="P")
                        {
                            std::vector<ID> chirAtms;
                            std::vector<ID> HAtms;
                            for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                                 iNB != iA->connAtoms.end(); iNB++)
                            {
                                
                                if (tAtoms[*iNB].chemType !="H")
                                {
                                    chirAtms.push_back(tAtoms[*iNB].id);
                                }
                                else
                                {
                                    HAtms.push_back(tAtoms[*iNB].id);
                                }
                            }
                            
                            int nH= (int)HAtms.size();
                            while((int)chirAtms.size() <3 && nH>0 )
                            {
                                int tPos = (int)HAtms.size() - nH; 
                                chirAtms.push_back(HAtms[tPos]);
                                nH--;
                            }
                        
                            if ((int)chirAtms.size() >=3)
                            {
                                std::string idxStr=IntToStr(idxC);
                                if (idxC <10)
                                {
                                    idxStr = "chir_0" + idxStr;
                                }
                                else
                                {
                                    idxStr = "chir_" + idxStr;
                                }
                                
                                // Not let H in as possible
                                  
                                outRestrF << longName 
                                          << std::setw(10) << idxStr 
                                          << std::setw(10)  << iA->id;
                       
                                outRestrF << std::setw(10) << chirAtms[0] << "    ";
                                outRestrF << std::setw(10) << chirAtms[1] << "    ";
                                outRestrF << std::setw(10) << chirAtms[2] << "    ";
                                outRestrF << std::setw(12) << "BOTH" << std::endl;
                                idxC++;
                            }
                        }
                    }
                }
                */
            }
            
            // Planar group section
            if ((int)tPlas.size() >0)
            {
                outRestrF << "loop_" << std::endl
                          << "_chem_comp_plane_atom.comp_id"  << std::endl
                          << "_chem_comp_plane_atom.plane_id" << std::endl
                          << "_chem_comp_plane_atom.atom_id"  << std::endl
                          << "_chem_comp_plane_atom.dist_esd" << std::endl;
                int idxP = 1;
                for (std::vector<PlaneDict>::iterator iP=tPlas.begin();
                        iP !=tPlas.end(); iP++)
                {
                    if (iP->atoms.size() > 3)
                    {
                        
                        std::string idxPStr = IntToStr(idxP);
                        idxPStr = "plan-" + idxPStr;
                        for(std::map<ID, int>::iterator iAt=iP->atoms.begin();
                                iAt != iP->atoms.end(); iAt++)
                        {
                            std::string tID;
                            if (iAt->first.find("\'") !=std::string::npos)
                            {
                                tID = "\"" + iAt->first + "\"";
                            }
                            else
                            {
                                tID = iAt->first;
                            }
                            outRestrF << longName
                                      << std::setw(10) << idxPStr
                                      << std::setw(12)  << tID
                                      << std::setw(8)  << "0.020" << std::endl;
                        }
                        idxP++;
                    }
                }
            }
            outRestrF.close();
        }    
    }
    
    
    extern void outMMCif2(FileName tFName, 
                         ID tMonoRootName,
                         ChemComp  &         tPropComp,
                         std::vector<LIBMOL::AtomDict>& tAtoms,
                         // std::vector<int>    & tHydroAtoms,
                         std::vector<LIBMOL::BondDict>& tBonds, 
                         std::vector<LIBMOL::AngleDict>& tAngs, 
                         std::vector<LIBMOL::TorsionDict>& tTorsions, 
                         std::vector<LIBMOL::RingDict> & tRings, 
                         std::vector<LIBMOL::PlaneDict>& tPlas, 
                         std::vector<LIBMOL::ChiralDict>& tChs)
    {
        
        for (std::vector<AtomDict>::iterator iA=tAtoms.begin();
                iA !=tAtoms.end(); iA++)
        {
            // std::cout << iA->id << std::endl;
            if (iA->id.find("\'") !=std::string::npos)
            {
                iA->id = "\"" + iA->id + "\"";
            }
            // std::cout << iA->id << std::endl;
        }
        
        // std::cout << "Print Pos " << std::endl;
        // Open a temp file for writing 
        std::vector<std::string> aSetStrs;
        StrTokenize(tFName, aSetStrs, '.');
        std::string outTempFName;
        for (unsigned i=0; i < aSetStrs.size()-1; i++ )
        {
            outTempFName+=aSetStrs[i];
        }
        outTempFName +="_ac.txt";
        std::cout << "output AandC file name : " << outTempFName << std::endl;
        
        std::ofstream outRestrF(tFName);
        
        if(outRestrF.is_open())
        {
            
            srand((unsigned)std::time( NULL ));
            // Temp 
            // 1. Global section 
            outRestrF << "global_" << std::endl
                    << "_lib_name         ?" << std::endl
                    << "_lib_version      ?" << std::endl
                    << "_lib_update       ?" << std::endl;
        
            
            
            // 'LIST OF MONOMERS' section
            
            
            std::string longName =tMonoRootName.substr(0,3);
            std::string sName =tMonoRootName.substr(0,3);
            
            //StrUpper(longName);
            
            
            ID ligType = "non-polymer";
            
            for (std::vector<LIBMOL::RingDict>::iterator iR=tRings.begin();
                    iR != tRings.end(); iR++)
            {   
                if (iR->sugarType.compare("pyranose")==0)
                {
                    ligType = "pyranose";
                    break;
                }
             
            }
            
            std::vector<ID>  aAATab;
            initAminoAcidTab(aAATab);
            if (isAminoAcid(aAATab, longName) && longName.find("PRO")==std::string::npos)
            {
                ligType = "L-peptide";
            }
            
   
            int nH = getHAtomNum(tAtoms);
            
            
            outRestrF << "# ------------------------------------------------" << std::endl
                    << "#" << std::endl
                    << "# ---   LIST OF MONOMERS ---" << std::endl
                    << "#" << std::endl
                    << "data_comp_list" << std::endl
                    << "loop_" << std::endl
                    << "_chem_comp.id" << std::endl
                    << "_chem_comp.three_letter_code" << std::endl
                    << "_chem_comp.name" << std::endl
                    << "_chem_comp.group" << std::endl
                    << "_chem_comp.number_atoms_all" << std::endl
                    << "_chem_comp.number_atoms_nh"  << std::endl
                    << "_chem_comp.desc_level" << std::endl;
            
            //if (tPropComp.id !=NullString)
            //{
            //    outRestrF << tPropComp.id  <<"\t"<< tPropComp.code << "\t" 
            //              << tPropComp.name << "\t" << tPropComp.group << "\t" 
            //              << tPropComp.numAtoms << "\t" 
            //              << tPropComp.numH << "\t."
                          // << (int)tAtoms.size()-(int)tHydroAtoms.size() << "\t."
            //              << std::endl;
                
            //}
            //else
            //{
            
            
            outRestrF << longName <<"\t"<< sName << "\t" << "'.\t\t'\t"
                      << ligType << "\t" << (int)tAtoms.size() << "\t" 
                      << (int)tAtoms.size()- nH << "\t."
                          // << (int)tAtoms.size()-(int)tHydroAtoms.size() << "\t."
                      << std::endl;
            
            outRestrF <<"# ------------------------------------------------------" << std::endl
                      <<"# ------------------------------------------------------" << std::endl
                      <<"#" << std::endl
                      <<"# --- DESCRIPTION OF MONOMERS ---" << std::endl
                      <<"#" << std::endl
                      <<"data_comp_" << longName << std::endl
                      <<"#" << std::endl; 
        
            if (tAtoms.size() >0)
            {
                // atom info section           
                outRestrF << "loop_" << std::endl
                          << "_chem_comp_atom.comp_id" << std::endl
                          << "_chem_comp_atom.atom_id" << std::endl
                          << "_chem_comp_atom.type_symbol" << std::endl
                          << "_chem_comp_atom.type_energy" << std::endl
                          << "_chem_comp_atom.charge" << std::endl
                          << "_chem_comp_atom.x" << std::endl
                          << "_chem_comp_atom.y" << std::endl
                          << "_chem_comp_atom.z" << std::endl;
                
                
                
                for (std::vector<AtomDict>::iterator iA = tAtoms.begin();
                        iA != tAtoms.end(); iA++)
                {
                    //double r1 =  (double) rand()/RAND_MAX;
                    //double r2 =  (double) rand()/RAND_MAX;
                    //double r3 =  (double) rand()/RAND_MAX;
                    REAL tCharge =0.0;
                    if (iA->charge !=0.0)
                    {
                        tCharge = iA->charge;
                    }
                    else if (iA->formalCharge !=0.0)
                    {
                        tCharge = iA->formalCharge;
                    }
                    
                    std::string strCharge = TrimSpaces(RealToStr(tCharge));
                    if (strCharge.find(".") !=strCharge.npos)
                    {
                        std::vector<std::string> tVec;
                        StrTokenize(strCharge, tVec, '.');
                        strCharge = tVec[0];
                    }
                    
                    StrUpper(iA->chemType);
                    StrUpper(iA->ccp4Type);
                    
                    outRestrF << longName
                              << std::setw(12) << iA->id 
                              << std::setw(6) << iA->chemType 
                              << std::setw(6) << iA->ccp4Type 
                              << std::setw(8) << strCharge 
                              << std::setw(12) << std::setprecision(3) << std::fixed 
                              << iA->coords[0] 
                              << std::setw(12) << std::setprecision(3) << std::fixed 
                              << iA->coords[1]
                              << std::setw(12) << std::setprecision(3) << std::fixed 
                              << iA->coords[2] << std::endl;
                    
                }
            }
            
            if (tBonds.size() >0)
            {
                // newly added 
                
                
                // Bond sections 
                outRestrF << "loop_" << std::endl
                          << "_chem_comp_bond.comp_id" << std::endl
                          << "_chem_comp_bond.atom_id_1" << std::endl
                          << "_chem_comp_bond.atom_id_2" << std::endl
                          << "_chem_comp_bond.type" << std::endl
                          << "_chem_comp_bond.aromatic"  << std::endl
                          << "_chem_comp_bond.value_dist"<< std::endl
                          << "_chem_comp_bond.value_dist_esd" << std::endl;
                          // << "_chem_comp_bond.exact_cod_dist" << std::endl;
               
                
                for (std::vector<BondDict>::iterator iB=tBonds.begin();
                          iB !=tBonds.end(); iB++)
                {
                    std::string tAr;
                    
                    
                    StrLower(iB->order);
                    
                    if (iB->order.find("arom") !=iB->order.npos)
                    {
                        tAr       = "y";
                    }
                    else
                    {
                        unifyStrForOrder(iB->order);
                        tAr = "n";
                    }
                                       
                    outRestrF <<  longName
                              << std::setw(12)  << tAtoms[iB->atomsIdx[0]].id  
                              << std::setw(12)  << tAtoms[iB->atomsIdx[1]].id  
                              << std::setw(12)  << iB->order 
                              << std::setw(8)   << tAr
                              << std::setw(10)  << std::setprecision(3)
                              << iB->value 
                              << std::setw(8) << std::setprecision(3)
                              << iB->sigValue << std::endl;
                    
                    //   if(iB->hasCodValue)
                    //   {
                    //       outRestrF << "Yes " << std::endl;
                    //   }
                    //   else
                    //   {
                    //       outRestrF << "No "  << std::endl;
                    //   }
                }
            }
           
            if (tAngs.size() > 0)
            {
                // Angle section
                outRestrF << "loop_" << std::endl
                          << "_chem_comp_angle.comp_id"   << std::endl
                          << "_chem_comp_angle.atom_id_1" << std::endl
                          << "_chem_comp_angle.atom_id_2" << std::endl
                          << "_chem_comp_angle.atom_id_3" << std::endl
                          << "_chem_comp_angle.value_angle"     << std::endl
                          << "_chem_comp_angle.value_angle_esd" << std::endl;
                          //  << "_chem_comp_angle.exact_cod_dist"  << std::endl;
                
                for (std::vector<AngleDict>::iterator iA=tAngs.begin();
                          iA != tAngs.end(); iA++)
                {
                    //for (std::vector<int>::iterator iAt=iA->atoms.begin();
                    //        iAt !=iA->atoms.end(); iAt++)
                    //{
                    // difference in comp atom definitions between cod and
                    // dictionary: inner-out1-out2(cod),
                    // atom1-atom2(center)-atom3(dictionary)
                    if (iA->sigValue < 1.5)
                    {
                        iA->sigValue = 1.5;
                    }
                    if (tAtoms[iA->atoms[0]].isMetal)
                    {
                        for (std::vector<REAL>::iterator iCA=iA->codAngleValues.begin();
                                iCA !=iA->codAngleValues.end(); iCA++)
                        {
                            outRestrF << tMonoRootName.substr(0,3) << std::setw(12)
                                      << tAtoms[iA->atoms[1]].id << std::setw(12)
                                      << tAtoms[iA->atoms[0]].id << std::setw(12)
                                      << tAtoms[iA->atoms[2]].id << std::setw(12) << std::setprecision(3) <<  *iCA << "    "
                                      << std::setw(8) << std::setprecision(2) << iA->sigValue << std::endl;
                        }
                    }
                    else
                    {
                        outRestrF << longName 
                                  << std::setw(12) << tAtoms[iA->atoms[1]].id
                                  << std::setw(12) << tAtoms[iA->atoms[0]].id 
                                  << std::setw(12) << tAtoms[iA->atoms[2]].id;
                        outRestrF << std::setw(12) << std::setprecision(3) <<  iA->value
                                  << std::setw(8) << std::setprecision(2) << iA->sigValue 
                                  << std::endl;
                    }
                    /*
                    if(iA->hasCodValue)
                     {
                       outRestrF << "Yes " << std::endl;
                     }
                     else
                     {
                       outRestrF << "No "  << std::endl;
                     }
                     */
                }
            }
        
            // Torsion section 
            if((int)tTorsions.size() !=0)
            {
                outRestrF << "loop_" << std::endl
                          << "_chem_comp_tor.comp_id"         << std::endl
                          << "_chem_comp_tor.id"              << std::endl
                          << "_chem_comp_tor.atom_id_1"       << std::endl
                          << "_chem_comp_tor.atom_id_2"       << std::endl
                          << "_chem_comp_tor.atom_id_3"       << std::endl
                          << "_chem_comp_tor.atom_id_4"       << std::endl
                          << "_chem_comp_tor.value_angle"     << std::endl
                          << "_chem_comp_tor.value_angle_esd" << std::endl
                          << "_chem_comp_tor.period"          << std::endl;
            
                int idxTor = 1;
                //std::cout << "number of torsions " << (int)allTorsions.size() << std::endl;
            
                for (std::vector<TorsionDict>::iterator iT=tTorsions.begin();
                        iT !=tTorsions.end(); iT++)
                {
                    //std::string idxTorStr=IntToStr(idxTor);
                    //idxTorStr = "tor_" + idxTorStr;
                    // std::cout << "Torsion angle " << idxTor 
                    //          << " It contains " << (int)iT->atoms.size() << std::endl;
                          
                    //std::cout << iT->atoms[0] << std::endl
                    //          << iT->atoms[1] << std::endl
                    //          << iT->atoms[2] << std::endl
                    //          << iT->atoms[3] << std::endl;
                
                    outRestrF << longName 
                              << std::setw(22) << iT->id
                              << std::setw(12)  << tAtoms[iT->atoms[0]].id 
                              << std::setw(12)  << tAtoms[iT->atoms[1]].id 
                              << std::setw(12)  << tAtoms[iT->atoms[2]].id 
                              << std::setw(12)  << tAtoms[iT->atoms[3]].id 
                              << std::setw(12) << std::setprecision(3) << iT->value  
                              << std::setw(8)  << "10.00" 
                              << std::setw(6)  << iT->period << std::endl;
                    idxTor++;        
                }
                
            }
            
            //  For chiral centers
            //bool l_ch = false;
            //if ((int)tChs.size() !=0)
            //{
            //    l_ch = true;
            //}
            //else
            //{
            //    for (int i_ch =0; i_ch < (int)tAtoms.size(); i_ch++)
            //    {
            //        if (tAtoms[i_ch].chiralIdx  > 0)
            //        {
            //            l_ch = true;
            //            break;
            //        }
            //    }
            //}
            
            //if (l_ch)
            if(tChs.size() !=0)
            {
                outRestrF << "loop_" << std::endl
                          << "_chem_comp_chir.comp_id" << std::endl
                          << "_chem_comp_chir.id" << std::endl
                          << "_chem_comp_chir.atom_id_centre" << std::endl
                          << "_chem_comp_chir.atom_id_1" << std::endl
                          << "_chem_comp_chir.atom_id_2" << std::endl
                          << "_chem_comp_chir.atom_id_3" << std::endl
                          << "_chem_comp_chir.volume_sign" << std::endl;
                
                // First the input chirals
                std::vector<ID>   inputChiralID;
                for (std::vector<ChiralDict>::iterator iCh = tChs.begin();
                        iCh != tChs.end(); iCh++)
                {
                   
                        inputChiralID.push_back(iCh->archID);
                        outRestrF << longName << "    " 
                                  << iCh->id  << "    ";
                        int numCh=0;
                        for (std::vector<int>::iterator iAt=iCh->atoms.begin();
                               iAt != iCh->atoms.end(); iAt++)
                        {
                            if (numCh < 4)
                            {
                                outRestrF << tAtoms[*iAt].id << "    ";
                                numCh++;
                            }
                        }
                        outRestrF << iCh->sign << std::endl;
                   
                }
                // New chiral that are not in the input list 
                /*
                int idxC =(int)inputChiralID.size();
                for (std::vector<AtomDict>::iterator iA=tAtoms.begin();
                          iA !=tAtoms.end(); iA++)
                {
                 
                    if (iA->chiralIdx == 1)
                    {
                        std::vector<ID>::iterator tFind;
                        tFind = std::find(inputChiralID.begin(), inputChiralID.end(), iA->id); 
                        if (tFind ==inputChiralID.end() && iA->chemType !="P")
                        {
                            std::vector<ID> chirAtms;
                            std::vector<ID> HAtms;
                            for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                                 iNB != iA->connAtoms.end(); iNB++)
                            {
                                
                                if (tAtoms[*iNB].chemType !="H")
                                {
                                    chirAtms.push_back(tAtoms[*iNB].id);
                                }
                                else
                                {
                                    HAtms.push_back(tAtoms[*iNB].id);
                                }
                            }
                            
                            int nH= (int)HAtms.size();
                            while((int)chirAtms.size() <3 && nH>0 )
                            {
                                int tPos = (int)HAtms.size() - nH; 
                                chirAtms.push_back(HAtms[tPos]);
                                nH--;
                            }
                        
                            if ((int)chirAtms.size() >=3)
                            {
                                std::string idxStr=IntToStr(idxC);
                                if (idxC <10)
                                {
                                    idxStr = "chir_0" + idxStr;
                                }
                                else
                                {
                                    idxStr = "chir_" + idxStr;
                                }
                                
                                // Not let H in as possible
                                  
                                outRestrF << longName 
                                          << std::setw(10) << idxStr 
                                          << std::setw(10)  << iA->id;
                       
                                outRestrF << std::setw(10) << chirAtms[0] << "    ";
                                outRestrF << std::setw(10) << chirAtms[1] << "    ";
                                outRestrF << std::setw(10) << chirAtms[2] << "    ";
                                outRestrF << std::setw(12) << "BOTH" << std::endl;
                                idxC++;
                            }
                        }
                    }
                }
                */
            }
            
            // Planar group section
            if ((int)tPlas.size() >0)
            {
                outRestrF << "loop_" << std::endl
                          << "_chem_comp_plane_atom.comp_id"  << std::endl
                          << "_chem_comp_plane_atom.plane_id" << std::endl
                          << "_chem_comp_plane_atom.atom_id"  << std::endl
                          << "_chem_comp_plane_atom.dist_esd" << std::endl;
                int idxP = 1;
                for (std::vector<PlaneDict>::iterator iP=tPlas.begin();
                        iP !=tPlas.end(); iP++)
                {
                    if (iP->atoms.size() > 3)
                    {
                        
                        std::string idxPStr = IntToStr(idxP);
                        idxPStr = "plan-" + idxPStr;
                        for(std::map<ID, int>::iterator iAt=iP->atoms.begin();
                                iAt != iP->atoms.end(); iAt++)
                        {
                            std::string tID;
                            if (iAt->first.find("\'") !=std::string::npos)
                            {
                                tID = "\"" + iAt->first + "\"";
                            }
                            else
                            {
                                tID = iAt->first;
                            }
                            outRestrF << longName
                                      << std::setw(10) << idxPStr
                                      << std::setw(12)  << tID
                                      << std::setw(8)  << "0.020" << std::endl;
                        }
                        idxP++;
                    }
                }
            }
            outRestrF.close();
        }    
    }
    
    extern void outMMCifFromOneMol(FileName tFName, 
                                   ID tMonoRootName,
                                   Molecule & tMol)
    {
        std::ofstream outRestrF(tFName);
        
        if(outRestrF.is_open())
        {
            ID rName;
            if (tMonoRootName.size() > 3)
            {
                rName = tMonoRootName.substr(0,3);
            }
            else
            {
                rName = tMonoRootName;
            }
            
            // 1. Global section 
            outRestrF << "global_" << std::endl
                    << "_lib_name         ?" << std::endl
                    << "_lib_version      ?" << std::endl
                    << "_lib_update       ?" << std::endl;
            
            if (tMol.atoms.size() !=0)
            {
                outRestrF << "loop_" << std::endl
                          << "_chem_comp_atom.comp_id" << std::endl 
                          << "_chem_comp_atom.atom_id" << std::endl  
                          << "_chem_comp_atom.type_symbol" << std::endl 
                          << "_chem_comp_atom.charge" << std::endl    
                          << "_chem_comp_atom.model_Cartn_x" << std::endl 
                          << "_chem_comp_atom.model_Cartn_y" << std::endl 
                          << "_chem_comp_atom.model_Cartn_z" << std::endl;
                
                for (std::vector<AtomDict>::iterator iAt=tMol.atoms.begin();
                        iAt !=tMol.atoms.end(); iAt++)
                {
                    StrUpper(iAt->chemType);
                    StrUpper(iAt->ccp4Type);
                    outRestrF << std::setw(8)<< rName
                              << std::setw(6) << iAt->id 
                              << std::setw(6) << iAt->chemType 
                              << std::setw(6) << iAt->ccp4Type 
                              << std::setw(8) << iAt->formalCharge 
                              << std::setw(12) << std::setprecision(3) << std::fixed 
                              << iAt->coords[0] 
                              << std::setw(12) << std::setprecision(3) << std::fixed 
                              << iAt->coords[1]
                              << std::setw(12) << std::setprecision(3) << std::fixed 
                              << iAt->coords[2] << std::endl;
                }
                
                std::cout << "#" << std::endl;
                
                if (tMol.bonds.size() !=0)
                {
                    outRestrF << "loop_" << std::endl
                              << "_chem_comp_bond.comp_id" << std::endl 
                              << "_chem_comp_bond.atom_id_1" << std::endl 
                              << "_chem_comp_bond.atom_id_2" << std::endl 
                              << "_chem_comp_bond.value_order" << std::endl
                              << "_chem_comp_bond.value_dist" << std::endl;
                    
                    for (std::vector<BondDict>::iterator iBo=tMol.bonds.begin();
                            iBo != tMol.bonds.end(); iBo++)
                    {
                        outRestrF << std::setw(8)<< rName
                                  << std::setw(6) << iBo->atoms[0] 
                                  << std::setw(6) << iBo->atoms[1] 
                                  << std::setw(6) << iBo->order 
                                  << std::setw(12) << std::setprecision(6) << std::fixed 
                                  << iBo->value << std::endl;
                    }
                }
                
                std::cout << "# ------------------------------------------------------" 
                          << std::endl;
                
            }
        }
    }
    
    extern int getHAtomNum(std::vector<LIBMOL::AtomDict>& tAtoms)
    {
        int tNumH=0;
        
        for (std::vector<AtomDict>::iterator iAt=tAtoms.begin();
                iAt !=tAtoms.end(); iAt++)
        {
            if (iAt->chemType.compare("H") ==0)
            {
                tNumH++;
            }
        }
        
        return tNumH;
    }
    
    extern void outHAtomIds(FileName tFName)
    {
        std::ofstream outHF(tFName);
    }
    
    extern void outMMCif3Secs(FileName tFName, 
                              ID tMonoRootName,
                              std::vector<LIBMOL::AtomDict>& tAtoms,
                              std::map<std::string, std::vector<std::string> > & tUnChangedEntries)
    {
       
        for (std::vector<AtomDict>::iterator iA=tAtoms.begin();
                iA !=tAtoms.end(); iA++)
        {
            if (iA->id.find("\'") !=std::string::npos)
            {
                iA->id = "\"" + iA->id + "\"";
            }
        }
        
        std::ofstream outRestrF(tFName);
        
        if(outRestrF.is_open())
        {
            if (tUnChangedEntries.find("head") !=tUnChangedEntries.end())
            {
                for (std::vector<std::string>::iterator iStr=tUnChangedEntries["head"].begin();
                        iStr != tUnChangedEntries["head"].end(); iStr++)
                {
                    outRestrF << *iStr;
                }
            }
            
            if (tAtoms.size() >0)
            {
                // atom info section
                std::string sName =tMonoRootName.substr(0,3);
                outRestrF << "loop_" << std::endl
                          << "_chem_comp_atom.comp_id" << std::endl
                          << "_chem_comp_atom.atom_id" << std::endl
                          << "_chem_comp_atom.type_symbol" << std::endl
                          << "_chem_comp_atom.type_energy" << std::endl
                          << "_chem_comp_atom.partial_charge" << std::endl
                          << "_chem_comp_atom.x" << std::endl
                          << "_chem_comp_atom.y" << std::endl
                          << "_chem_comp_atom.z" << std::endl;
                
                for (std::vector<AtomDict>::iterator iA = tAtoms.begin();
                        iA != tAtoms.end(); iA++)
                {
                    //double r1 =  (double) rand()/RAND_MAX;
                    //double r2 =  (double) rand()/RAND_MAX;
                    //double r3 =  (double) rand()/RAND_MAX;
                    StrUpper(iA->chemType);
                    StrUpper(iA->ccp4Type);
                    
                    std::string tId;
                    if (iA->id.find("\'") !=std::string::npos)
                    {
                      tId = "\"" + iA->id + "\"";
                    }
                    else
                    {
                        tId = iA->id;
                    }
           
                    // std::cout << "here " << tId << std::endl;
                    outRestrF << sName 
                              << std::setw(12) << tId 
                              << std::setw(6) << iA->chemType 
                              << std::setw(6) << iA->ccp4Type 
                              << std::setw(8) << iA->parCharge 
                              << std::setw(12) << std::setprecision(3) << std::fixed 
                              << iA->coords[0] 
                              << std::setw(12) << std::setprecision(3) << std::fixed 
                              << iA->coords[1]
                              << std::setw(12) << std::setprecision(3) << std::fixed 
                              << iA->coords[2] << std::endl;
                }
            }
            
            if (tUnChangedEntries.find("others") !=tUnChangedEntries.end())
            {
                for (std::vector<std::string>::iterator iStr=tUnChangedEntries["others"].begin();
                        iStr != tUnChangedEntries["others"].end(); iStr++)
                {
                    outRestrF << *iStr;
                }
            }
        }
    }
    
    extern void outAtomTypesAndConnections(FileName tFName,
                                        std::vector<LIBMOL::AtomDict>& tAtoms,
                                        std::vector<LIBMOL::BondDict>& tBonds,
                                        std::vector<LIBMOL::RingDict> & tRings)
    {
        
        if (tAtoms.size() !=0 && tBonds.size())
        {
            std::ofstream outTempF(tFName);
            if (outTempF.is_open())
            {
                outTempF << "ATOMS:" << std::endl;
                int nA=0;
                std::map<int, int> seriMap;
                std::vector<int> hAtomIdxs;
                for (std::vector<AtomDict>::iterator iA = tAtoms.begin();
                        iA != tAtoms.end(); iA++)
                {
                    outTempF << std::setw(10) << nA 
                             << std::setw(8)  << iA->chemType
                             << std::setw(10) << iA->id 
                             << std::setw(iA->codClass.size()+6) 
                             << iA->codClass << std::endl;
                    seriMap[iA->seriNum] = nA;
                    
                    if (iA->chemType.compare("H")==0)
                    {
                        hAtomIdxs.push_back(nA);
                    }
                    
                    nA++;
                    
                    
                }
                
                if (hAtomIdxs.size() >0)
                {
                    outTempF << "H FormType:" << std::endl;
                    for (std::vector<int>::iterator iH = hAtomIdxs.begin();
                        iH != hAtomIdxs.end(); iH++)
                    {
                        if (tAtoms[*iH].formType.size()==2 
                            && tAtoms[*iH].connAtoms.size() ==1)
                        {
                            int cIdx = tAtoms[*iH].connAtoms[0];
                            std::string cId = tAtoms[cIdx].id;
                            outTempF << std::setw(10) << tAtoms[*iH].formType[1]
                                     << std::setw(8)  << tAtoms[*iH].id
                                     << std::setw(8) <<  cId << std::endl;
                        }
                    }
                }
                outTempF << "CONNECTIONS:" << std::endl;
                for (std::vector<BondDict>::iterator iB=tBonds.begin();
                          iB !=tBonds.end(); iB++)
                {
                    outTempF << std::setw(12) << seriMap[tAtoms[iB->atomsIdx[0]].seriNum] 
                             << std::setw(12) << seriMap[tAtoms[iB->atomsIdx[1]].seriNum]
                             << std::endl;
                }
                
                if (tRings.size() > 0)
                {
                    outTempF << "Ring Information: " << std::endl;
                    int  idxR = 1;
                    std::string aR = "RING_";
                    for (std::vector<RingDict>::iterator iR=tRings.begin();
                          iR !=tRings.end(); iR++)
                    {
                        std::string aro ="";
                        if (iR->isAromatic)
                        {
                            aro = "Aromatic";
                        }
                        else if(iR->isAntiAroma)
                        {
                            aro = "Anti-Aromatic";
                        }
                        else
                        {
                            aro = "Non-Aromatic";
                        }
                        
                        std::string aLab = aR + IntToStr(idxR);
                        idxR++;
                        for (std::vector<AtomDict>::iterator 
                                iA=iR->atoms.begin();
                                iA !=iR->atoms.end(); iA++)
                        {
                            outTempF << std::setw(12) << aLab
                                     << std::setw(12) << iA->id
                                     << std::setw(20) << aro
                                     << std::endl;
                        }
                        
                    }
                }
                
                outTempF.close();  
            } 
        }
    }
    
    extern void outMetalAtomInfo(FileName tFName,
                                 GenCifFile  & tCifObj)
    {
        std::ofstream outMF(tFName);
        
        outMF << "Is the structure not from Power diff : " 
              << tCifObj.notPowder << std::endl;
        outMF << "Resolution OK : " 
              << tCifObj.resolOK << std::endl;
        outMF << "R factors OK : " 
              << tCifObj.RFactorOK << std::endl;
        outMF << "Occp OK : " << tCifObj.colidOK << std::endl;
        outMF << "The following are metal atoms in the file " << std::endl;
        for (std::vector<AtomDict>::iterator iAt=tCifObj.allAtoms.begin();
                iAt != tCifObj.allAtoms.end(); iAt++)
        {
            if (iAt->isMetal)
            {
                outMF << "Atom : " << iAt->id << " : " 
                      << iAt->chemType << std::endl;
            }
        }
        
        outMF.close();        
              
    }
    
    void extern outSelectedAtomInfo(FileName tFName,
                                    std::vector<std::string> & tSelectedIds)
    {
        std::ofstream outAF(tFName);
        
        
        outAF << "The following are selected atoms in the file " << std::endl;
        for (std::vector<std::string>::iterator iI = tSelectedIds.begin();
             iI != tSelectedIds.end(); iI++)
        {              
            outAF << "Atom : " << *iI << std::endl;    
        }
        
        outAF.close();
    }
    
}
