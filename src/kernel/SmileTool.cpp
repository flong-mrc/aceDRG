/* 
 * File:   SmileTool.cpp
 * Author: flong
 *
 * Created on September 11, 2012, 1:55 PM
*/

#include "SmileTool.h"

namespace LIBMOL
{
    SmileTool::SmileTool():itsCurStr(NullString),
                           itsPrevStr(NullString), 
                           itsCurPos(ZeroInt),
                           itsPrevPos(ZeroInt),
                           itsCurAtomIdx(ZeroInt),
                           itsCurAtom(NullPoint),
                           itsPrevAtom(NullPoint)
    {
        setOrgSet();
        setSpecialChars();
    }
    
    SmileTool::SmileTool(std::string& tSStr):itsCurStr(NullString),
                           itsPrevStr(NullString), 
                           itsCurPos(ZeroInt),
                           itsPrevPos(ZeroInt),
                           itsCurAtomIdx(ZeroInt),
                           itsCurAtom(NullPoint),
                           itsPrevAtom(NullPoint)
    {
         setOrgSet();
         setSpecialChars();
         SmileToMol(tSStr);
    }
    
    
    
    SmileTool::SmileTool(const std::vector<AtomDict>& tAllAtoms, 
                         const std::vector<BondDict>& tAllBonds, 
                         const std::map<ID,std::vector<RingDict> >& tAllRings)
    :itsCurStr(NullString),itsPrevStr(NullString),
     itsCurPos(ZeroInt), itsPrevPos(ZeroInt),
     itsCurAtomIdx(ZeroInt),itsCurAtom(NullPoint),
     itsPrevAtom(NullPoint)
    {
        setOrgSet();
        setSpecialChars();
    }
    
    SmileTool::~SmileTool()
    {
        if(itsCurAtom)
        {
            delete itsCurAtom;
            itsCurAtom = NULL;
        }
        if(itsPrevAtom)
        {
            delete itsPrevAtom;
            itsPrevAtom = NULL;
        }
    }
    
    void SmileTool::setOrgSet()
    {
        OrganSet["C"].push_back(4);
        OrganSet["N"].push_back(3);
        OrganSet["N"].push_back(5);
        OrganSet["O"].push_back(2);
        OrganSet["S"].push_back(2);
        OrganSet["S"].push_back(4);
        OrganSet["S"].push_back(6);
        OrganSet["P"].push_back(3);
        OrganSet["P"].push_back(5);
        OrganSet["B"].push_back(3);
        OrganSet["F"].push_back(1);
        OrganSet["Cl"].push_back(1);
        OrganSet["Br"].push_back(1);
        OrganSet["I"].push_back(1);  
    }
    
    void SmileTool::setSpecialChars()
    {
        // atoms
        specialChars["atom"].push_back("[");
        specialChars["atom"].push_back("]");
        specialChars["atom"].push_back("+");
        specialChars["atom"].push_back("-");
        
        // Bonds
        specialChars["bond"].push_back("-");    // or just ""
        specialChars["bond"].push_back("=");
        specialChars["bond"].push_back("#");
        specialChars["bond"].push_back(":");
        
        
        // Branch 
        specialChars["branch"].push_back("(");
        specialChars["branch"].push_back(")");
        specialChars["branch"].push_back("%");
        
        // Disconnected Structures
        specialChars["discStru"].push_back(".");
        
        // chiral center
        specialChars["chiral"].push_back("\\");
        specialChars["chiral"].push_back("/");
        specialChars["chiral"].push_back("@");
        specialChars["chiral"].push_back("@@");
         
    }
    
    bool SmileTool::isOrgAtom(ID tId)
    {
        bool aT = false;
        
        if (OrganSet.find(tId) !=OrganSet.end())
        {
            aT = true;
        }
        
        return aT;
    }
    
    void SmileTool::SmileToMol(ID & tSS, Molecule & tMol)
    {
        
        if ((int)tSS.size() !=0)
        {
            itsCurPos =0;
            
            while (itsCurPos < (int)tSS.size())
            {
                ID id1=tSS.substr(itsCurPos, 1);
                if (id1.compare("[")==0)
                {
                    AtomDict aAtom;
                    itsCurPos++;
                    
                    //while(id1.compare("]")!=0 && itsCurPos < (int)tSS.size())
                    //{
                    //    id1 = tSS.substr(itsCurPos, 1);
                    //    if(isInt(id1) )
                    //    {
                            
                    //    }
                    //    itsCurPos++;
                    //}
                    //setOneAtomProp(tSS, itsCurPos);
                }
                if (std::find(specialChars["atom"].begin(),
                              specialChars["atom"].end(), id1) !=specialChars["atom"].end())
                {
                    setOneAtomProp(tSS, id1);
                }
                else if (std::find(specialChars["bond"].begin(),
                              specialChars["bond"].end(), id1)!=specialChars["atom"].end())
                {
                    setOneBondProp(tSS, id1);
                }
                else if  (std::find(specialChars["branch"].begin(),
                              specialChars["branch"].end(), id1)!=specialChars["atom"].end())
                {
                    setBranch(tSS, id1);
                }
                else if (std::find(specialChars["chiral"].begin(),
                              specialChars["chiral"].end(), id1)!=specialChars["atom"].end())
                {
                    setChiral(tSS, id1);
                }
                else if (specialChars["discStru"][0].compare(id1)==0)
                {
                    setDisconnectedUnits(tSS, id1);
                }
                else
                {       
 
                    // 
                }
                
            }
        }
    }
    
    void SmileTool::SmileToMol(ID & tSS)
    {
        
    }
    
    void SmileTool::checkAndSetOneAtom(ID tSS, ID tId)
    {
        /*
        if (isValidAtomName(tId))
        {
            if (itsCurPos + 1 < (int)tSS.size())
            {
                ID id2=tSS.substr(itsCurPos+1, 1);
                ID id3=tId+id2;
                if (isValidAtomName(id3))
                {
                    setOneAtomName(tSS, id3);
                }
                            else
                            {
                                setOneAtomName(tSS, tId);
                            }
                        }
                        else
                        {
                            setOneAtomName(tSS, tId);
                        }
                    }
                    else
                    {
                        ID id4 = tSS.substr(itsCurPos, 1);
                        StrUpper(id4);
                        
                        if (isValidAtomName(id4))
                        {
                            if (itsCurPos + 1 < (int)tSS.size())
                            {
                                ID id2=tSS.substr(itsCurPos+1, 1);
                                ID id3=id4+id2;
                                
                                if (isValidAtomName(id3))
                                {
                                    setOneAtomName(tSS, id3);
                                }
                                else
                                {
                                    setOneAtomName(tSS, id4);
                                }
                            }
                            else
                            {
                                setOneAtomName(tSS, id4);
                            }
                        }
                  }
                         
         }
    */
                           
    }
   
    void SmileTool::setOneAtomName(ID & tSS, ID tId)
    {
        AtomDict tAtom;
        ID id1 = tId.substr(0,1);
        StrUpper(id1);
        if (id1.compare(tId.substr(0,1)) ==0)
        {
            tAtom.chemType=tId;
            itsPrevPos = itsCurPos;
            itsCurPos++;
        }
        else
        {
            if((int)tId.size() ==2)
            {
                tAtom.chemType=id1+tId.substr(1,1);
            }
            else
            {
                tAtom.chemType=id1;
            }
            // the atom is in  an aromatic ring
            tAtom.bondingIdx = 2;
            tAtom.chiralIdx  = -1;
            itsPrevPos = itsCurPos;
            itsCurPos+=2;
        }
            
        tAtom.cChemType  = tId;
                
        // allAtoms.push_back(tAtom);
        itsCurAtomIdx++;
    }
    
    void SmileTool::setOneAtomProp(ID& tSS, ID tId)
    {
    }
    void SmileTool::setBranch(ID& tSS, ID tId)
    {
    }
    
    void SmileTool::setChiral(ID& tSS, ID tId)
    {
    }
    
    void SmileTool::setDisconnectedUnits(ID& tSS, ID tId)
    {
    }
    
    void SmileTool::setOneBondProp(ID & tSS, ID tId)
    {
        
            
            if (tId.compare("=")==0)
            {
                BondDict aBond;
                itsPrevPos = itsCurPos;
                itsCurPos+=1;
                setOneAtomName(tSS, tSS.substr(itsCurPos, 2));
                
                
                // aBond.fullAtoms[allAtoms[itsCurAtomIdx-1].chemType] = itsCurAtomIdx-1; 
            }
    }
    
    void SmileTool::SmileToCifUsingLibcheck(FileName tSmiFileName,
                                            std::string & tOutFileName )
    {
        std::string clibcheck(std::getenv("CBIN"));
        if (clibcheck.empty())
        {
            std::cerr << "You need setup CCP4 suite first " << std::endl;
            exit(1);
        }
        clibcheck.append("/libcheck");
        
        std::string tSmiStrs(tSmiFileName);
        std::string tLogName(tSmiFileName);
        tLogName.append("_libcheck.log");
        std::string tOutName(tSmiFileName);
        tOutName.append("_libcheck_out");
        std::string aCmdLine(clibcheck);
        aCmdLine.append("  << eof > " + tLogName + " \n\n");
        aCmdLine.append("FILE_SMILE  ");
        aCmdLine.append(tSmiStrs + "\n");
        aCmdLine.append("FILE_O  ");
        aCmdLine.append(tOutName + "\n");
        aCmdLine.append("eof \n");
        aCmdLine.append("\n");
        std::cout << "Comline is " << aCmdLine << std::endl;
        
        cmdExecute(aCmdLine);
        
        tOutFileName.clear();
        
        tOutFileName.append(tOutName + ".lib");
        
        // should be outside the function
        if (!isFileExist(tOutFileName.c_str()))
        {
            std::cout << "libcheck failed to generate lib file "
                    << tOutFileName << std::endl;
            tOutFileName.clear();
            exit(1);
        }   
    }
     
    
}