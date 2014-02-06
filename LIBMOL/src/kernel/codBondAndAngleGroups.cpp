/* 
 * File:   codBondAndAngleGroups.h
 * Author: flong
 *
 * Created on August 30, 2012, 5:31 PM
 */

#include "codBondAndAngleGroups.h"

namespace LIBMOL
{
    CodBondAndAngleGroup::CodBondAndAngleGroup()
    {
    }
    
    CodBondAndAngleGroup::CodBondAndAngleGroup(FileName tBFName=NULL, 
                              FileName tAFName=NULL)
    {
        
        groupCodBonds(tBFName);
        
        groupCodAngles(tAFName);
        
    }
    
    CodBondAndAngleGroup::~CodBondAndAngleGroup()
    {
    }
    
    void CodBondAndAngleGroup::groupCodBonds(FileName tBFName) 
                                       
    {
        try
        {
            std::cout << "Clustering all COD bonds " << std::endl;
            
            
            std::ifstream codBondFile(tBFName);
            if(codBondFile.is_open())
            {
                std::string tRecord="";
                int nline = 0;
                while(!codBondFile.eof())
                {
                    std::getline(codBondFile, tRecord);
                    tRecord = TrimSpaces(tRecord);
                    std::vector<std::string> tBuf;
                    StrTokenize(tRecord, tBuf);
                    
                    if ((int)tBuf.size() ==11)
                    {
                     
                        int ha1, ha2;
                        
                        ha1 = StrToInt(TrimSpaces(tBuf[0]));
                        ha2 = StrToInt(TrimSpaces(tBuf[1]));
                        
                        
                        for (int i=2; i < (int)tBuf.size(); i++)
                        {
                            tBuf[i] = TrimSpaces(tBuf[i]);
                        }
                       
                        allDictBondsIdx[ha1][ha2][tBuf[2]][tBuf[3]][tBuf[4]][tBuf[5]][tBuf[6]][tBuf[7]]= nline;
                        
                        
                        BondDict aBond;
                        aBond.seriNum = nline;
                        aBond.atomsHashingCodes.push_back(ha1);
                        aBond.atomsHashingCodes.push_back(ha2);
                        aBond.atomsNBRep.push_back(tBuf[4]);
                        aBond.atomsNBRep.push_back(tBuf[5]);
                        aBond.atomsCodClasses.push_back(tBuf[6]);
                        aBond.atomsCodClasses.push_back(tBuf[7]);
                        aBond.value    = StrToReal(tBuf[8]);
                        aBond.sigValue = StrToReal(tBuf[9]);
                        if(aBond.sigValue < 0.01)
                        {
                            aBond.sigValue = 0.01;
                        }
                        aBond.numCodValues = StrToInt(tBuf[10]);
                        
                        allDictBonds.push_back(aBond);
                        
                        nline+=1;
                        
                    }
                }
                
                codBondFile.close();
            }
        }
        catch (std::exception & e)
        {
            std::cout << e.what() << std::endl;
        }
        std::cout << "Finish grouping COD bonds " << std::endl;
       
    }
    
    void CodBondAndAngleGroup::groupCodAngles(FileName tAFName)
    {
        try
        {
            // should be something like std::string tNewCodBondFileName(clibMonDir + "/list/bonds.txt");
            std::ifstream codAngleFile(tAFName);
            int nline =0;
            if(codAngleFile.is_open())
            {
                std::string tRecord="";
                int iLine = 0;
        
                while(!codAngleFile.eof())
                {
                    std::getline(codAngleFile, tRecord);
                    tRecord = TrimSpaces(tRecord);
                    std::vector<std::string> tBuf;
                    StrTokenize(tRecord, tBuf);
                    
                    if ((int)tBuf.size() ==15)
                    {
                        int ha1, ha2, ha3;
                        
                        ha1 = StrToInt(TrimSpaces(tBuf[0]));
                        ha2 = StrToInt(TrimSpaces(tBuf[1]));
                        ha3 = StrToInt(TrimSpaces(tBuf[2]));
                        
                        for (int i=3; i < (int)tBuf.size(); i++)
                        {
                            tBuf[i] = TrimSpaces(tBuf[i]);
                        }
                        
                        allDictAnglesIdx[ha1][ha2][ha3][tBuf[3]][tBuf[4]][tBuf[5]][tBuf[6]][tBuf[7]][tBuf[8]][tBuf[9]][tBuf[10]][tBuf[11]]
                                =nline;

                        
                        AngleDict  aAngle;
                        aAngle.seriNum    = iLine;
                        iLine++;
                        aAngle.atomsNBRep.push_back(tBuf[6]);
                        aAngle.atomsNBRep.push_back(tBuf[7]);
                        aAngle.atomsNBRep.push_back(tBuf[8]);
                        
                        aAngle.atomsCodClasses.push_back(tBuf[9]);
                        aAngle.atomsCodClasses.push_back(tBuf[10]);
                        aAngle.atomsCodClasses.push_back(tBuf[11]);
                        
                        aAngle.value        = StrToReal(tBuf[12]);
                        aAngle.sigValue     = StrToReal(tBuf[13]);
                        if(aAngle.sigValue <0.0001)
                        {
                            aAngle.sigValue = 3.0;
                        }
                        aAngle.numCodValues = StrToInt(tBuf[14]);
                        
                        allDictAngles.push_back(aAngle);
                        
                        nline +=1;
                    }
                }
                
                codAngleFile.close();
            
            }
        }
        catch (std::exception & e)
        {
            std::cout << e.what() << std::endl;
        }
        
    }
}