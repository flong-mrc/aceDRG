/* 
 * File:   ExtraRestrDictFile.h
 * Author: flong
 *
 * Created on September 22, 2011, 5:42 PM
 */

#include "ExtraRestrDictFile.h"

namespace LIBMOL
{
    // -----Member functions in RestraitListSet -------
   
    RestraitListSet::RestraitListSet():label(NullString),
            groupLabel(NullString)
    {
    }
 
    RestraitListSet::RestraitListSet(const RestraitListSet & tP)
    {
        label      = tP.label;
        groupLabel = tP.groupLabel;
 
        if(tP.monoSet.size())
        {
            for (std::vector<RestrOneMonoInfo>::const_iterator iT
                    =tP.monoSet.begin();
                    iT != tP.monoSet.end(); iT++)
            {
                monoSet.push_back(*iT);
            }
        }
        
        if (tP.restrBondList.size()) 
        {
            for (std::vector<Bond>::const_iterator iT=tP.restrBondList.begin();
                    iT != tP.restrBondList.end(); iT++)
            {
                restrBondList.push_back((*iT));
            }
        }

        if (tP.restrAngleList.size())
        {          
            for (std::vector<Angle>::const_iterator iT=tP.restrAngleList.begin();
                    iT != tP.restrAngleList.end(); iT++)
            {
                restrAngleList.push_back((*iT));
            }
        }
        
        if (tP.restrTorsionList.size())
        {
            for (std::vector<Torsion>::const_iterator iT=tP.restrTorsionList.begin();
                    iT != tP.restrTorsionList.end(); iT++)
            {
                restrTorsionList.push_back((*iT));
            }
        }
        
        if (tP.restrChiralList.size())
        {
            for (std::vector<Chiral>::const_iterator iT= tP.restrChiralList.begin();
                    iT != tP.restrChiralList.end(); iT++)
            {
                restrChiralList.push_back((*iT));
            }
        }
        
        if (tP.restrPlaneList.size())
        {
            for (std::vector<Plane>::const_iterator iT=tP.restrPlaneList.begin();
                    iT != tP.restrPlaneList.end(); iT++)
            {
                restrPlaneList.push_back((*iT));
            }
        }
    }
    
    RestraitListSet::~RestraitListSet()
    {
        if (! monoSet.empty())
        {
            monoSet.clear();
        }
        if (!restrBondList.empty())
        {
            restrBondList.clear();
        }
        
        if (!restrAngleList.empty())
        {
            restrAngleList.clear();
        }   
        
        if (!restrTorsionList.empty())
        {
            restrTorsionList.clear();
        }
        
        if (!restrChiralList.empty())
        {
            restrChiralList.clear();
        }
        
        if (!restrPlaneList.empty())
        {
            restrPlaneList.clear();
        }        
    }
    
    // -----Member functions in ExtraRestrDictFile -------
    
    ExtraRestrDictFile::ExtraRestrDictFile ()   
    {
        
        itsName               = NullString;
        itsFullPath           = NullString; 
        itsExt                = NullString;
        itsSize               = ZeroInt;
        itsOpenMode           = std::ios::in;
        itsType               = UNKNOWN; 
        itsState.isOpen       = false;
        itsState.isReadable   = true; 
        itsState.isWritable   = false; 
        itsState.isExecutable = false;
        itsInfoType           = DICT_ADDTIONAL;
        itsCurOneMono         = NULL;
        itsCurRestrListSet       = NULL;
    }
    
    ExtraRestrDictFile::ExtraRestrDictFile (Name  tFName)
    {
        itsName               = NullString;
        itsFullPath           = NullString; 
        itsExt                = NullString;
        itsSize               = ZeroInt;
        itsOpenMode           = std::ios::in;
        itsType               = UNKNOWN; 
        itsState.isOpen       = false;
        itsState.isReadable   = true; 
        itsState.isWritable   = false; 
        itsState.isExecutable = false;
        itsInfoType           = DICT_ADDTIONAL;
        itsCurOneMono         = NULL;
        itsCurRestrListSet       = NULL;
       
        inFile.open(tFName.c_str(), std::ios::in);
        setupSystem();
    }
    
    ExtraRestrDictFile::ExtraRestrDictFile(FileName tFName)
    {
        itsName               = tFName;
        itsFullPath           = NullString; 
        itsExt                = NullString;
        itsSize               = ZeroInt;
        itsOpenMode           = std::ios::in;
        itsType               = UNKNOWN; 
        itsState.isOpen       = false;
        itsState.isReadable   = true; 
        itsState.isWritable   = false; 
        itsState.isExecutable = false;
        itsInfoType           = DICT_ADDTIONAL;
        itsCurOneMono         = NULL;
        itsCurRestrListSet       = NULL;
        inFile.open(tFName, std::ios::in);
   
        setupSystem();
       
               
      
        /*
        for (int i=0; i < allRestrLists.restrBondList.size(); i++)
        {
            std::cout << "Bond "  << i << "  : "
                      <<allRestrLists.restrBondList[i].getLength(true) 
                      <<std::endl;
           std::cout << "sigma    : "
                      <<allRestrLists.restrBondList[i].getSigLength() 
                      <<std::endl;
        }
        
        
        std::cout << "Number of  torsion restraints : " 
                << allRestrLists.restrTorsionList.size() << std::endl;
       
        for (int i=0; i < allRestrLists.restrTorsionList.size(); i++)
        {
            std::cout << "Torsion "  << i << "  : "
                      <<allRestrLists.restrTorsionList[i].getValue(true) 
                      <<std::endl;
            std::cout << "sigma    : "
                      <<allRestrLists.restrTorsionList[i].getSigValue() 
                      <<std::endl;
            std::cout << "Numb of Atoms    : "
                      <<allRestrLists.restrTorsionList[i].atoms.size() 
                      <<std::endl;
        }      
       
        
        std::cout << "Number of chiral restraints : " 
                << allRestrLists.restrChiralList.size() << std::endl;
        
        for (int i=0; i < allRestrLists.restrChiralList.size(); i++)
        {
            std::cout << "Chiral "  << i << "  : "
                      <<allRestrLists.restrChiralList[i].getValue(true) 
                      <<std::endl;
            std::cout << "sigma    : "
                      <<allRestrLists.restrChiralList[i].getSigValue() 
                      <<std::endl;
            std::cout << "Numb of Atoms    : "
                      <<allRestrLists.restrChiralList[i].atoms.size() 
                      <<std::endl;
        }       
        */
        
    }
    
    ExtraRestrDictFile::~ExtraRestrDictFile ()
    {
        if(itsCurRestrListSet !=NULL)
        {
            delete itsCurRestrListSet;
            itsCurRestrListSet = NULL;
        }
        if(itsCurOneMono != NULL)
        {
            delete itsCurOneMono;
            itsCurOneMono = NULL;
        }
    }
    
    RestrInfoType ExtraRestrDictFile::getRecordType(ALine tLine)
    {
        int placeHolder=0;
        std::cout << "Just list placeHolder " << placeHolder << std::endl;
        return DICT_PAIRS;   
    }
    
    
    void ExtraRestrDictFile::setupSystem()
    {
    
        if (inFile.is_open())
        {    
            bool tOK = true;
            std::string tRecord="";
  
            while(!inFile.eof() && tOK)
            {
                std::getline(inFile, tRecord);
                if (inFile.good())
                {   
                    std::vector<std::string> tStrGrps;
                    StrTokenize(tRecord, tStrGrps);
                    // std::cout << tRecord << std::endl;
                    
                    if (tStrGrps.size() !=0)
                    {
                        std::string tHeader= TrimSpaces(tStrGrps[0]);
                        if(!tHeader.compare("addrestraints"))
                        {
                            itsInfoType = DICT_ADDTIONAL;
                            tOK = setAddRestrDefInfo(tRecord);
                            
                        }
                        if (!tHeader.compare("pairs"))
                        {
                            itsInfoType = DICT_PAIRS;
                            tOK = setPairRestrDefInfo(tRecord);
                        }
                        if (!tHeader.compare("torsion"))
                        {
                            tOK = setOneTorsionRestr(tRecord);
                        }
                        if (!tHeader.compare("bond"))
                        {
                            
                            tOK = setOneBondRestr(tRecord);
                        }
                        if (!tHeader.compare("angle"))
                        {
                            tOK = setOneAngleRestr(tRecord);
                        }
                        if (!tHeader.compare("chiral"))
                        {
                            tOK = setOneChiralRestr(tRecord);
                        }
                        if (!tHeader.compare("plane"))
                        {
                            tOK = setOnePlaneRestr(tRecord);
                        }
                    }
                }
            } 
            
           
            if(itsCurRestrListSet !=NULL)
            {
                allRestrLists.push_back(*itsCurRestrListSet);
                delete itsCurRestrListSet;
                itsCurRestrListSet = NULL;
            }
            
            inFile.close();
            
            // Check what the restraints from the dictionary file are like
            // std::cout<< "Finish read restraints from the dictionary"
            //        << std::endl;
            std::cout << std::endl
                      <<"------------------------------------------" 
                      << std::endl;
            
            if ((int)itsInfoType == 0)
            {
                std::cout << "Type of restraints to be found: " << std::endl
                          << "additional restraints " << std::endl;
            }
            else if ((int)itsInfoType == 1)
            {
                std::cout << "Type of restraints to be found: " << std::endl
                          << "base pair restraints " << std::endl;
            }
            else
            {
                std::cout << "Type of restraints to be found: " << std::endl
                          << "Unknown " << std::endl;
            }
            
            std::cout <<"------------------------------------------" 
                      << std::endl;
            
            if(allRestrLists.size() != 0)
            {
                
                for (int i=0; i<(int) allRestrLists.size(); i++  )
                {
                    std::cout << "Dictionary value for restraint set : "
                              << i+1 << std::endl;
                    for (int j=0; j < (int)allRestrLists[i].monoSet.size(); j++)
                    {
                        std::cout << "Monomer: " 
                                << allRestrLists[i].monoSet[j].monomerName
                                << "  Group: " 
                                << allRestrLists[i].monoSet[j].groupName 
                                << std::endl;
                    }
                    std::cout <<"Label: " << allRestrLists[i].label << std::endl;
                    std::cout << "Group Label: " << allRestrLists[i].groupLabel
                              << std::endl;
                    //std::cout << "Number of bond restraints " 
                    //          <<  allRestrLists[i].restrBondList.size()
                    //          << std::endl;
                    for (int j=0; j < (int)allRestrLists[i].restrBondList.size(); j++)
                    {
                        std::cout << "Bond "  << j+1 << "  : " << std::endl
                        << "Atom : " 
                        << allRestrLists[i].restrBondList[j].atoms[0].getName()
                        << " Atom : "
                        << allRestrLists[i].restrBondList[j].atoms[1].getName()
                        << " Value "
                        <<allRestrLists[i].restrBondList[j].getLength(true) 
                        << " sigma "
                        <<allRestrLists[i].restrBondList[j].getSigLength() 
                        <<std::endl;
                    }
                    
                    //std::cout << "Number of Torsions "
                    //        << allRestrLists[i].restrTorsionList.size()
                    //        << std::endl;
                    for (int j=0; j < (int)allRestrLists[i].restrTorsionList.size(); 
                            j++)
                    {
                        std::cout << "Torsion " << j+1 << " : " << std::endl
                        << "Atom "
                        << allRestrLists[i].restrTorsionList[j].atoms[0].getName()
                        << " Atom "
                        << allRestrLists[i].restrTorsionList[j].atoms[1].getName()
                        << " Atom "
                        << allRestrLists[i].restrTorsionList[j].atoms[2].getName()
                        << " Atom "
                        << allRestrLists[i].restrTorsionList[j].atoms[3].getName()
                        << " value " 
                        << allRestrLists[i].restrTorsionList[j].getValue(true)
                        << " sigma "
                        << allRestrLists[i].restrTorsionList[j].getSigValue()
                        << std::endl;
                        
                    }
                    
                    for (int j=0; j < (int)allRestrLists[i].restrChiralList.size(); 
                            j++)
                    {
                        std::cout << "Chiral " << j+1 << " : " << std::endl
                        << "Atom "
                        << allRestrLists[i].restrChiralList[j].atoms[0].getName()
                        << " Atom "
                        << allRestrLists[i].restrChiralList[j].atoms[1].getName()
                        << " Atom "
                        << allRestrLists[i].restrChiralList[j].atoms[2].getName()
                        << " Atom "
                        << allRestrLists[i].restrChiralList[j].atoms[3].getName()
                        << " value " 
                        << allRestrLists[i].restrChiralList[j].getValue(true)
                        << " sigma "
                        << allRestrLists[i].restrChiralList[j].getSigValue()
                        << std::endl;
                        
                    }
                    
                    std::cout <<"------------------------------------------" 
                              << std::endl;
                }
                std::cout << std::endl;
            }    
        }
        else
        {
            std::cerr << itsName << "could not be opened for reading "
                    << std::endl;
            exit(1);
        }
        
    }
    
    bool ExtraRestrDictFile::setAddRestrDefInfo(ALine tLine)
    {
        // At the moment, it is the same as setPairRestrDefInfo
        // May change it in the future
        return setPairRestrDefInfo(tLine);
        
    }
    
    bool ExtraRestrDictFile::setPairRestrDefInfo(ALine  tLine)
    {
        std::vector<std::string> curStrs;
        StrTokenize(tLine, curStrs);

        if (curStrs.size())
        {   
           
            if (itsCurRestrListSet == NULL)
            {
                itsCurRestrListSet = new RestraitListSet();
            }
            else
            {     
                allRestrLists.push_back(*itsCurRestrListSet);
                delete itsCurRestrListSet;
                itsCurRestrListSet = NULL;
                itsCurRestrListSet = new RestraitListSet();
            }
        
            for (int i = 0; i < (int)curStrs.size()-1; i++)
            {
                if (curStrs[i].find("monomer") != std::string::npos)
                {
                    itsCurOneMono = new RestrOneMonoInfo;
                    itsCurOneMono->monomerName = TrimSpaces(curStrs[i+1]);
               
                }
                if (curStrs[i].compare("group") ==0)
                {
                    itsCurOneMono->groupName = TrimSpaces(curStrs[i+1]);
                
                    itsCurRestrListSet->monoSet.push_back(*itsCurOneMono);
                  
                    delete itsCurOneMono;
                    itsCurOneMono = NULL;    
                }         
                if (curStrs[i].size() == 5 
                        && curStrs[i].find("label") != std::string::npos)
                {
                    itsCurRestrListSet->label = TrimSpaces(curStrs[i+1]);
                }               
                if (curStrs[i].find("grouplabel") != std::string::npos)
                {
                   
                    itsCurRestrListSet->groupLabel = TrimSpaces(curStrs[i+1]);            
                }
            }
            curStrs.clear();
            return true;
        }
        return false;
    }
 
    bool ExtraRestrDictFile::setOneTorsionRestr(ALine  tLine)
    {
        std::vector<std::string> curStrs;
        StrTokenize(tLine, curStrs);
        if (curStrs.size())
        {
            Torsion aTorsion;
            for (int i = 0; i < (int)curStrs.size()-1; i++)
            {
                if (curStrs[i].find("atom") != std::string::npos)
                {
                    Atom aAtom;
                    aAtom.setName(curStrs[i+1]);
                    aAtom.setResName(curStrs[i+2]);
                    aTorsion.atoms.push_back(aAtom);
                }
                if (curStrs[i].find("value") != std::string::npos)
                {
                    aTorsion.setValue(StrToReal(curStrs[i+1]), true);
                }              
                if (curStrs[i].find("sigma") != std::string::npos)
                {
                    aTorsion.setSigValue(StrToReal(curStrs[i+1]));
                }
            }
            itsCurRestrListSet->restrTorsionList.push_back(aTorsion);
            return true;
        }
        return false;
    }
    
    bool ExtraRestrDictFile::setOneBondRestr(ALine  tLine)
    {
        
        std::vector<std::string> curStrs;
        StrTokenize(tLine, curStrs);
        if (curStrs.size())
        {
            Bond aBond;
            
            for (int i = 0; i < (int)curStrs.size()-1; i++)
            {
                if (curStrs[i].find("atom") != std::string::npos)
                {
                    Atom aAtom;
                    aAtom.setName(TrimSpaces(curStrs[i+1]));
                    aAtom.setResName(curStrs[i+2]);
                    // std::cout << aAtom.getResName() << std::endl;
                    aBond.atoms.push_back(aAtom);
                }
                if (curStrs[i].find("value") != std::string::npos)
                {
                    aBond.setLength(StrToReal(curStrs[i+1]), true);
                }              
                if (curStrs[i].find("sigma") != std::string::npos)
                {
                    aBond.setSigLength(StrToReal(curStrs[i+1]));
                }
            }
            
            itsCurRestrListSet->restrBondList.push_back(aBond);
            
            /* for (int i =0; i < itsCurRestrListSet->restrBondList.size();
                    i++)
            {
                std::cout << "atom 1 ";
                std::cout << itsCurRestrListSet->restrBondList[i].atoms[0].getName()
                        << std::endl;
                std::cout << "atom 2 ";
                std::cout << itsCurRestrListSet->restrBondList[i].atoms[1].getName()
                        << std::endl;
                std::cout << "bond length ";
                std::cout << itsCurRestrListSet->restrBondList[i].getLength(true) 
                        << std::endl;
            }
             */
            
            return true;
        }
        return false;        
    }
    
    bool ExtraRestrDictFile::setOneAngleRestr(ALine  tLine)
    {
        std::vector<std::string> curStrs;
        StrTokenize(tLine,curStrs);
        if (curStrs.size())
        {
            Angle aAngle;
            for (int i = 0; i < (int)curStrs.size()-1; i++)
            {
                if (curStrs[i].find("atom") != std::string::npos)
                {
                    Atom aAtom;
                    aAtom.setName(curStrs[i+1]);
                    aAtom.setResName(curStrs[i+2]);
                    aAngle.atoms.push_back(aAtom);
                }
                if (curStrs[i].find("value") != std::string::npos)
                {
                    aAngle.setValue(StrToReal(curStrs[i+1]), true);
                }              
                if (curStrs[i].find("sigma") != std::string::npos)
                {
                    aAngle.setSigValue(StrToReal(curStrs[i+1]), true);
                }
            }
            itsCurRestrListSet->restrAngleList.push_back(aAngle);
            return true;
        }
        return false;            
    }
    
    bool ExtraRestrDictFile::setOneChiralRestr(ALine  tLine)
    {
        std::vector<std::string> curStrs;
        StrTokenize(tLine, curStrs);
        if (curStrs.size())
        {
            
            Chiral aChiral;
            for (int i = 0; i < (int)curStrs.size()-1; i++)
            {
                if (curStrs[i].find("atom") != std::string::npos)
                {
                    Atom aAtom;
                    aAtom.setName(curStrs[i+1]);
                    aAtom.setResName(curStrs[i+2]);
                    aChiral.atoms.push_back(aAtom);
                }
                if (curStrs[i].find("value") != std::string::npos)
                {
                    aChiral.setValue(StrToReal(curStrs[i+1]), true);
                }              
                if (curStrs[i].find("sigma") != std::string::npos)
                {
                    aChiral.setSigValue(StrToReal(curStrs[i+1]));
                }
            }
            itsCurRestrListSet->restrChiralList.push_back(aChiral);
            return true;
        }
        return false;            
    }
    
    bool ExtraRestrDictFile::setOnePlaneRestr(ALine  tLine)
    {
        std::vector<std::string> curStrs;
        StrTokenize(tLine, curStrs);
        if (curStrs.size())
        {
            Plane aPlane;
            for (int i = 0; i < (int)curStrs.size()-1; i++)
            {
                if (curStrs[i].find("atom") != std::string::npos)
                {
                    Atom aAtom;
                    aAtom.setName(curStrs[i+1]);
                    aAtom.setResName(curStrs[i+2]);
                    aPlane.atoms.push_back(aAtom);
                }
                if (curStrs[i].find("value") != std::string::npos)
                {
                    aPlane.setValue(StrToReal(curStrs[i+1]), true);
                }              
                if (curStrs[i].find("sigma") != std::string::npos)
                {
                    aPlane.setSigValue(StrToReal(curStrs[i+1]));
                }
            }
            itsCurRestrListSet->restrPlaneList.push_back(aPlane);
            return true;
        }
        return false;            
    }
    
    RestrInfoType ExtraRestrDictFile::getRestrInfoType() 
    {
        return itsInfoType;
    }
    void ExtraRestrDictFile::setRestrInfoType(RestrInfoType  tR)
    {
        itsInfoType = tR;
    }
    
 
    
    
    
}
