/* 
 * File:   spaceGroupTable.cpp
 * Author: flong
 *
 * Created on August 7, 2013, 6:46 PM
 */


#include "crystInfo.h"

namespace LIBMOL
{
    SpaceGroupMember::SpaceGroupMember():sgNum(-1)
    {
    }
    
    SpaceGroupMember::SpaceGroupMember(const SpaceGroupMember & tSG):sgNum(tSG.sgNum)
    {
        for(std::map<ID, std::vector<ID> >::const_iterator iS=tSG.sgSymb.begin();
                iS !=tSG.sgSymb.end(); iS++)
        {
            for (std::vector<ID>::const_iterator iSymb=iS->second.begin();
                    iSymb !=iS->second.end(); iSymb++)
            {
                sgSymb[iS->first].push_back(*iSymb);
            }
        }
        
        for (std::map<std::string, std::vector<std::vector<REAL> > >::const_iterator iO=tSG.sgOp.begin();
                iO != tSG.sgOp.end(); iO++)
        {
            sgOp[iO->first]=iO->second;
        }
    }
    
    SpaceGroupMember::~SpaceGroupMember()
    {
    }
    
    SpaceGroupTable::SpaceGroupTable()
    {   
    }
    
    SpaceGroupTable::~SpaceGroupTable()
    {   
        if (itsCurSG !=NULL)
        {
            delete itsCurSG;
            itsCurSG = NULL;
        }
    }
    
    void SpaceGroupTable::getSgtableFromCCP4()
    {
        std::string clibDir(std::getenv("CLIBD"));
        std::string symInfoName = clibDir + "/syminfo.lib";
        std::ifstream symInfoFile(symInfoName.c_str());
        int i=0;
        if(symInfoFile.is_open())
        {
            std::string tRecord="";
            while(!symInfoFile.eof())
            {
                std::getline(symInfoFile, tRecord);    
                tRecord = TrimSpaces(tRecord);
                if (tRecord.find("begin_spacegroup") !=std::string::npos)
                {
                    if (i==0)
                    {
                        i++;
                    }
                    else
                    {
                        allSgs.push_back(*itsCurSG);
                        delete itsCurSG;
                        itsCurSG = NULL;
                    }
                    
                    itsCurSG= new SpaceGroupMember();
                }
                else if (itsCurSG !=NULL)
                {
                    TrimSpaces(tRecord);
                    std::vector<std::string> tBuf;
                    
                    if (tRecord.find("number") !=std::string::npos) 
                    {
                        StrTokenize(tRecord, tBuf);
                        if ((int)tBuf.size()==2)
                        {
                            itsCurSG->sgNum = StrToInt(tBuf[1]);
                        }
                    }
                    else if (tRecord.find("symbol") !=std::string::npos)
                    {
                        if(tRecord.find("xHM") !=std::string::npos)
                        {
                            StrTokenize(tRecord, tBuf, '\'');
                            for (int j=1; j <(int)tBuf.size(); j++)
                            {
                                if ((int)tBuf[j].size()!=0)
                                {
                                    itsCurSG->sgSymb["xHM"].push_back(TrimSpaces(tBuf[j]));
                                }
                            }
                        }
                        else if(tRecord.find("Hall") !=std::string::npos)
                        {
                            StrTokenize(tRecord, tBuf, '\'');
                            for (int j=1; j <(int)tBuf.size(); j++)
                            {
                                if ((int)tBuf[j].size()!=0)
                                {
                                    itsCurSG->sgSymb["Hall"].push_back(TrimSpaces(tBuf[j]));
                                }
                            }
                        }
                        else if(tRecord.find("laue") !=std::string::npos)
                        {
                            StrTokenize(tRecord, tBuf, '\'');
                            for (int j=1; j <(int)tBuf.size(); j++)
                            {
                                if ((int)tBuf[j].size()!=0)
                                {
                                    itsCurSG->sgSymb["laue"].push_back(TrimSpaces(tBuf[j]));
                                }
                            }
                        }
                        else if(tRecord.find("symop") !=std::string::npos
                                || tRecord.find("cenop") !=std::string::npos )
                        {
                            StrTokenize(tRecord, tBuf);
                            
                            if ((int)tBuf.size() ==2 && 
                                 tBuf[1] !="x,y,z")
                            {
                                std::vector<std::string>          tBuf2;
                                StrTokenize(tBuf[1], tBuf2, ',');
                                std::vector<std::vector<REAL> >   aMat;
                                StrToSymmOps(tBuf2, aMat);
                                itsCurSG->sgOp[tBuf[1]]=aMat;
                            }
                        }
                    }
                }
            }
        }
    }
    
    Cell::Cell():a(ZeroReal),
            b(ZeroReal),
            c(ZeroReal),
            alpha(ZeroReal),
            beta(ZeroReal),
            gamma(ZeroReal),
            vol(ZeroReal),
            lattice(NullString)
    { 
    }
            
    Cell::Cell(const Cell & tC):a(tC.a),
            b(tC.b),
            c(tC.c),
            alpha(tC.alpha),
            beta(tC.beta),
            gamma(tC.gamma),
            vol(tC.vol),
            lattice(tC.lattice)
    { 
    }
    
    Cell::~Cell()
    {
    }
 
   
    Resolution::Resolution(): resolLimit(RESOLTHRESHOLD), 
                              dMax(-1.0),
                              thetaMax(0.0),
                              wavLen(-1.0),
                              lSet(false)

    {
    }
    
    Resolution::Resolution(const Resolution & tResol): resolLimit(tResol.resolLimit),
                                                       dMax(tResol.dMax),
                                                       thetaMax(tResol.thetaMax),
                                                       wavLen(tResol.wavLen),
                                                       lSet(tResol.lSet)
    {
    }
    
    Resolution::~Resolution()
    {
    }
    
    
    void Resolution::setResol()
    {
        if (wavLen >0.0 && fabs(thetaMax) > 1.0 )
        {
            dMax =wavLen/(2*sin(thetaMax*PI/180.0));
            lSet = true;
        }
    }
    
    void Resolution::setResolLimit(double tResolLimit)
    {
        resolLimit = tResolLimit;  
    }
    
    
    CrystInfo::CrystInfo() : itsCell(NullPoint),
                             itsSpaceGroup(NullPoint),
                             itsResolution(NullPoint)
    {           
    }
    
    CrystInfo::CrystInfo(const CrystInfo & tCryst)
    {
        itsCell = new Cell(*tCryst.itsCell);
        itsSpaceGroup = new SpaceGroupMember(*tCryst.itsSpaceGroup);
        itsResolution = new Resolution(*tCryst.itsResolution);
        
        for (std::vector<REAL>::const_iterator iO=tCryst.ORIGXn.begin();
                iO != tCryst.ORIGXn.end(); iO++)
        {
            ORIGXn.push_back(*iO);
        }
        for (std::vector<REAL>::const_iterator iM=tCryst.MTRIXn.begin();
                iM != tCryst.MTRIXn.end(); iM++)
        {
            MTRIXn.push_back(*iM);
        }
        
        for (std::vector<REAL>::const_iterator iS=tCryst.SCALEn.begin();
                iS != tCryst.SCALEn.end(); iS++)
        {
            SCALEn.push_back(*iS);
        }
        
    }
    
    CrystInfo:: ~CrystInfo()
    {
        if(itsCell)
        {
            delete itsCell;
            itsCell = NullPoint;
        }
            
        if(itsSpaceGroup)
        {
            delete itsSpaceGroup;
            itsSpaceGroup = NullPoint;
        }
        
        if(itsResolution)
        {
            delete itsResolution;
            itsResolution = NullPoint;
        }
    }
    
}

