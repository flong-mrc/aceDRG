/* 
 * File:   MonomerLib.h
 * Author: flong
 *
 * Created on September 27, 2011, 02:32 PM
 */

#include "MonomerLib.h"

namespace LIBMOL
{
    MonomerLib::MonomerLib():itsMonoListName(NullPoint)
    {
    }
    
    MonomerLib::MonomerLib(FileName  tF)
    {
        setupSystem(tF);
    }
    
    MonomerLib::~MonomerLib()
    {
    }
    
    FileName  MonomerLib::getMonoListName()
    {
        return itsMonoListName;
    }
    void MonomerLib::setMonoListName(FileName tF)
    {
        itsMonoListName = tF;
    }
    
    void MonomerLib::setupSystem(FileName  tF)
    {
        if(tF)
        {
            setMonoListName(tF);
            std::ifstream tDictFile(itsMonoListName);
            //std::cout << "Dictionary monomer List file : " << std::endl
            //        << itsMonoListName << std::endl;
            
            if (tDictFile.is_open())
            {
               
                std::string tRecord="";
            
                bool lTable = false;
                while( !tDictFile.eof())
                {
                    std::getline(tDictFile, tRecord);
                    
                    if (tDictFile.good())
                    { 
                        std::vector<std::string> tStrs; 
                        StrTokenize(tRecord, tStrs);
                       // std::cout << tRecord << std::endl;
                       
                        if (tRecord.find("LIST OF RECOGNIZED MODIFIED MONOMERS") 
                                != std::string::npos)
                        {
                            lTable = false;
                        }
                        if (tRecord.find("_chem_comp.desc_level") 
                                != std::string::npos)
                        {
                            lTable = true;
                        }  
                        if(lTable &&  tRecord.length() == 79 
                                && tStrs.size() > 2)
                        {
                            std::string tKey =  TrimSpaces(tStrs[1]);
                            std::vector<std::string> tStrs2;
                            StrTokenize(TrimSpaces(tRecord.substr(52)), tStrs2);
                            monomerGroups[tKey] = TrimSpaces(tStrs2[0]);
                        }
                    }
                }
                tDictFile.close();
                std::cout << "------------------------------------------"
                          << std::endl;
                std::cout<<"Read monomer Lib, its size is : " 
                        << (int)monomerGroups.size() << std::endl;
                std::cout << "------------------------------------------"
                          << std::endl;
                /*for (std::map<std::string, std::string>::const_iterator 
                        iT=monomerGroups.begin();
                        iT != monomerGroups.end(); iT++)
                {
                    std::cout << "Monomer Name: " << (*iT).first << std::endl;
                    std::cout << "Monomer Group: " << (*iT).second << std::endl;
                }
                 */
            }
            else 
            {
                throw setMonomerException();
            }
        }
    } 
    
}
    