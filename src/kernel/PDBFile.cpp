
/* 
 * File:   PDBFile.cpp
 * Author: flong
 *
 * Created on September 2, 2011, 11:07 AM
 */

#include "PDBFile.h"

#ifndef FILE_H
#include "file.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef LIBG_ATOMASSEMBLY_H
#include "atomAssembly.h"
#endif

#ifndef RESIDUE_H
#include "residue.h"
#endif

#ifndef CHAIN_H
#include "chain.h"
#endif

#ifndef LIBG_MODEL_H
#include "libgmodel.h"
#endif

#ifndef SECONDARYSTRUCTURES_H
#include "secondaryStructures.h"
#endif



namespace LIBMOL
{
    
    PDBFile::PDBFile() : itsTempRecordType(UnknownType),
                         itsCurrentModel(NullPoint),
                         itsCurrentAtom(NullPoint),
                         itsCurrentAltLoc(false),
                         itsCurrentHetAtm(NullPoint),
                         itsCurrentResidue(NullPoint),
                         itsCurrentModRes(NullPoint),
                         itsCurrentChain(NullPoint),
                         itsCurrentHelix(NullPoint),
                         itsCurrentSheet(NullPoint),
                         itsCurrentSSBond(NullPoint),
                         itsCrystInfo(NullPoint)
    {       
        allHeaderInfo.CAVEAT= NULL;
        allHeaderInfo.COMPND= NULL;
        allHeaderInfo.HEADER= NULL;
        allHeaderInfo.KEYWDS= NULL;
        allHeaderInfo.OBSLTE= NULL;
        allHeaderInfo.SOURCE= NULL;
        allHeaderInfo.SPLIT= NULL;
        allHeaderInfo.TITLE= NULL;
    }
    
    PDBFile::PDBFile(Name               tFname,
            std::ios_base::openmode tOpenMode) :
            itsTempRecordType(UnknownType),
            itsCurrentModel(NullPoint),
            itsCurrentAtom(NullPoint),
            itsCurrentAltLoc(false),
            itsCurrentHetAtm(NullPoint),
            itsCurrentResidue(NullPoint),
            itsCurrentModRes(NullPoint),
            itsCurrentChain(NullPoint),
            itsCurrentHelix(NullPoint),
            itsCurrentSheet(NullPoint),
            itsCurrentSSBond(NullPoint),
            itsCrystInfo(NullPoint)
    {
        if (tOpenMode == std::ios::in)
        {
            inFile.open(tFname.c_str(), tOpenMode);
            itsCurrentModel   = new Model();
            itsCurrentAtom    = new Atom();
            itsCurrentAltLoc  = false;
            itsCurrentHetAtm  = new Atom();
            itsCurrentResidue = new Residue();
            itsCurrentModRes  = new ModRes();
            itsCurrentChain   = new Chain();
            itsCurrentHelix   = new Helix();
            itsCurrentSheet   = new Sheet();
            itsCurrentSSBond  = new SSBond();
            itsCrystInfo      = new CrystInfo();
            setupSystem();
        }
        else
        {
            outFile.open(tFname.c_str(), tOpenMode);
        }

    }
    
   PDBFile::PDBFile(FileName        tFname,
            std::ios::openmode      tOpenMode=std::ios::in ) :
            itsTempRecordType(UnknownType),
            itsCurrentModel(NullPoint),
            itsCurrentAtom(NullPoint),
            itsCurrentAltLoc(false),
            itsCurrentHetAtm(NullPoint),
            itsCurrentResidue(NullPoint),
            itsCurrentModRes(NullPoint),
            itsCurrentChain(NullPoint),
            itsCurrentHelix(NullPoint),
            itsCurrentSheet(NullPoint),
            itsCurrentSSBond(NullPoint),
            itsCrystInfo(NullPoint)
    {
        if (tOpenMode == std::ios::in)
        {       
            inFile.open(tFname, tOpenMode);
            itsCurrentModel   = new Model();
            itsCurrentAtom    = new Atom();
            itsCurrentHetAtm  = new Atom();
            itsCurrentResidue = new Residue();
            itsCurrentModRes  = new ModRes();
            itsCurrentChain   = new Chain();
            itsCurrentHelix   = new Helix();
            itsCurrentSheet   = new Sheet();
            itsCurrentSSBond  = new SSBond();
            itsCrystInfo      = new CrystInfo();            
            setupSystem();
        
            // all check
            
            //std::cout << "There are " << allAtomList.size() 
            //        << "  atoms " << std::endl;
           
            //std::cout<< "There are(is) " << allModels.size() 
            //        << " model in the system" << std::endl;
            //for (int i= 0; i < (int)allModels.size(); i++)
            //{
                //std::cout << "in model: " << i 
                //        << ", there are " << allModels[i].chains.size() 
                //        << " chains " << std::endl;
                //for (int j =0; j < (int)allModels[i].chains.size(); j++ )
                //{
                    //std::cout<< "for Chain " << j+1
                    //        << " : its ID " << allModels[i].chains[j].getID() 
                    //        << std::endl;
                    //std::cout << "There are " << allModels[i].chains[j].residues.size()
                    //        << " residues" << std::endl;
                    
                    //for (int k =0; k < (int)allModels[i].chains[j].residues.size();
                    //        k++)
                    // {   
                        
                        //std::cout << "for Residue " 
                        //        << allModels[i].chains[j].residues[k].getName()
                        //       << "  :  " << std::endl;
                        //std::cout << "Its chain ID " 
                        //        << allModels[i].chains[j].residues[k].getChainID()
                        //        << std::endl;
                        //std::cout << "Its seqNum " 
                        //        << allModels[i].chains[j].residues[k].getSeqNum()
                        //        << std::endl;
                        
                        //std::cout << "There are " 
                        //        <<  allModels[i].chains[j].residues[k].atoms.size()
                        //        << " atoms. " << std::endl;
                        //for (int l = 0; l < (int)allModels[i].chains[j].residues[k].atoms.size();
                        //        l++)
                        //{
                        //    if(allModels[i].chains[j].residues[k].atoms[l].altLoc.size()!=0)
                        //    {
                        //    std::cout << " atom " 
                        //            << allModels[i].chains[j].residues[k].atoms[l].getSeriNum()
                        //            << "  :  " << std::endl;
                        //    std::cout << " Name " 
                        //            << allModels[i].chains[j].residues[k].atoms[l].getName()
                        //            << std::endl;
                        //    std::cout << "Its seqNum " 
                        //              << allModels[i].chains[j].residues[k].atoms[l].getSeqNum()
                        //              << std::endl;
                        //    std::cout << " Chain ID " 
                        //            << allModels[i].chains[j].residues[k].atoms[l].getChainID()
                        //            << std::endl;
                        //    std::cout << "Coordinate X : "
                        //            << allModels[i].chains[j].residues[k].atoms[l].coords[0]
                        //            << std::endl;
                        //    std::cout << "AltLocation Coordinate X : "
                        //            << allModels[i].chains[j].residues[k].atoms[l].altLoc[0]
                        //            << std::endl;
                        //    }
                        
                        //}
                          
                    //}
                //}
            //}
        }
        else
        {
            outFile.open(tFname, tOpenMode);
        }
    }
         
    PDBFile::~PDBFile()
    {
        if(inFile.is_open())
        {
            inFile.close();
        }
        
        if(outFile.is_open())
        {
            outFile.close();
        }
        // delete operations
    }
    
    void PDBFile::setPDBRecordHeader()
    {
        PDBRecordHeader.push_back("UNKNOWN");
        PDBRecordHeader.push_back("ANISOU");
        PDBRecordHeader.push_back("ATOM  ");
        PDBRecordHeader.push_back("AUTHOR");
        PDBRecordHeader.push_back("CAVEAT");
        PDBRecordHeader.push_back("CISPEP");
        PDBRecordHeader.push_back("COMPND");
        PDBRecordHeader.push_back("CONECT");
        PDBRecordHeader.push_back("CRYST1");
        PDBRecordHeader.push_back("DBREF ");
        PDBRecordHeader.push_back("DBREF1");
        PDBRecordHeader.push_back("DBREF2");
        PDBRecordHeader.push_back("END   ");
        PDBRecordHeader.push_back("ENDMDL");
        PDBRecordHeader.push_back("EXPDTA");
        PDBRecordHeader.push_back("FORMUL");
        PDBRecordHeader.push_back("FTNOTE");
        PDBRecordHeader.push_back("HEADER");
        PDBRecordHeader.push_back("HELIX ");
        PDBRecordHeader.push_back("HET   ");
        PDBRecordHeader.push_back("HETATM");
        PDBRecordHeader.push_back("HETNAM");
        PDBRecordHeader.push_back("HETSYN");
        PDBRecordHeader.push_back("HYDBND");
        PDBRecordHeader.push_back("JRNL  ");
        PDBRecordHeader.push_back("KEYWDS");
        PDBRecordHeader.push_back("LINK  ");
        PDBRecordHeader.push_back("MASTER");
        PDBRecordHeader.push_back("MODEL ");
        PDBRecordHeader.push_back("MODRES");
        PDBRecordHeader.push_back("MTRIX1");
        PDBRecordHeader.push_back("MTRIX2");
        PDBRecordHeader.push_back("MTRIX3");
        PDBRecordHeader.push_back("OBSLTE");
        PDBRecordHeader.push_back("ORIGX1");
        PDBRecordHeader.push_back("ORIGX2");
        PDBRecordHeader.push_back("ORIGX3");
        PDBRecordHeader.push_back("REMARK");
        PDBRecordHeader.push_back("REVDAT");
        PDBRecordHeader.push_back("SCALE1");
        PDBRecordHeader.push_back("SCALE2");
        PDBRecordHeader.push_back("SCALE3");
        PDBRecordHeader.push_back("SEQADV");
        PDBRecordHeader.push_back("SEQRES");
        PDBRecordHeader.push_back("SHEET ");
        PDBRecordHeader.push_back("SIGATM");
        PDBRecordHeader.push_back("SIGUIJ");
        PDBRecordHeader.push_back("SITE  ");
        PDBRecordHeader.push_back("SLTBRG");
        PDBRecordHeader.push_back("SOURCE");
        PDBRecordHeader.push_back("SPRSDE");
        PDBRecordHeader.push_back("SSBOND");
        PDBRecordHeader.push_back("TER   ");
        PDBRecordHeader.push_back("TITLE ");
        PDBRecordHeader.push_back("TURN  ");
        PDBRecordHeader.push_back("TVECT ");
    }
       /*   "DBREF1",
            "DBREF2",
            "END   ",
            "ENDMDL",
            "EXPDTA",
            "FORMUL",
            "FTNOTE",
            "HEADER",
            "HELIX ",
            "HET   ",
            "HETATM",
            "HETNAM",
            "HETSYN",
            "HYDBND",
            "JRNL  ",
            "KEYWDS",
            "LINK  ",
            "MASTER",
            "MODEL ",
            "MODRES",
            "MTRIX1",
            "MTRIX2",
            "MTRIX3",
            "OBSLTE",
            "ORIGX1",
            "ORIGX2",
            "ORIGX3",
            "REMARK",
            "REVDAT",
            "SCALE1",
            "SCALE2",
            "SCALE3",
            "SEQADV",
            "SEQRES",
            "SHEET ",
            "SIGATM",
            "SIGUIJ",
            "SITE  ",
            "SLTBRG",
            "SOURCE",
            "SPRSDE",
            "SSBOND",
            "TER   ",
            "TITLE ",
            "TURN  ",
            "TVECT "); */
           
    void PDBFile::setupSystem()
    {
          
        if (inFile.is_open() )
        {    
            setPDBRecordHeader();

            bool tOK = true;
            std::string tRecord="";
            
            while(!inFile.eof() && tOK)
            {
                std::getline(inFile, tRecord);
               
                if (inFile.good())
                {   
                    int tRecType = getRecordType(tRecord); 
                    
                    switch(tRecType)
                    {
                        /*
                        case HeadType:
                        {
                            tOK=extractTitleInfo(tRecord);
                            break;
                        }
                        case SeqResType:  // Primary Structure Section
                        {
                            tOK = extractSEQRESInfo(tRecord);
                            break;
                        }
                        case ModResType:
                        {
                            tOK = extractMODRESInfo(tRecord);
                            break;
                        }
                        case HelixType:   // Secondary Structure Section
                        {
                            tOK = extractHelixInfo(tRecord);
                            break;
                        }
                        case SheetType:
                        {
                            tOK = extractSheetInfo(tRecord);
                            break;
                        }
                        case SsBondType:  // Connectivity Annotation Section
                        {
                            tOK=extractSSBondInfo(tRecord);
                            break;
                        }
                        case LinkType:
                        {
                            tOK=extractLinkInfo(tRecord);
                            break;
                        }
                        case Cryst1Type:   // Crystallographic Section
                        {
                            tOK = extractCryst1Info(tRecord);
                            break;
                        }
                        case ORIGX1Type:
                        {
                            tOK = extractOrigXNInfo(tRecord);
                            break;
                        }
                        case ORIGX2Type:
                        {
                            tOK = extractOrigXNInfo(tRecord);
                            break;
                        }
                        case ORIGX3Type:
                        {
                            tOK = extractOrigXNInfo(tRecord);
                            break;
                        }
                        case MTRIX1Type:
                        {
                            tOK = extractMatrixNInfo(tRecord);
                            break;
                        }
                        case MTRIX2Type:
                        {
                            tOK = extractMatrixNInfo(tRecord);
                            break;
                        }
                        case MTRIX3Type:
                        {
                            tOK = extractMatrixNInfo(tRecord);
                            break;
                        }
                        case SCALE1Type:
                        {
                            tOK = extractScaleNInfo(tRecord);
                            break;
                        }
                        case SCALE2Type:
                        {
                            tOK = extractScaleNInfo(tRecord);
                            break;
                        }
                        case SCALE3Type:
                        {
                            tOK = extractScaleNInfo(tRecord);
                            break;
                        }
                        */
                        case ModelType:     // coordinate section
                        {
                            initOneModel(tRecord);
                            break;
                        }
                        case AtomType:
                        {
                            tOK=extractAtomInfo(tRecord);
                            break;
                        }
                        case AnisouType:
                        {
                            tOK=extractAnisouInfo(tRecord);
                            break;
                        }
                        case TerType:
                        {
                            tOK=extractTerInfo(tRecord);
                            break;
                        }
                        /*
                        case HetatmType:
                        {
                            tOK=extractAtomInfo(tRecord);
                            break;
                        }
                         */
                        case EndMdlType:
                        {
                            tOK=extractEndMdlInfo(tRecord);
                            break;
                        }
                        /*
                        case ConectType:
                        {
                            tOK=extractConnectInfo(tRecord);
                            break;
                        }
                       */
                        default :
                        {
                            itsTempRecordType=UnknownType;
                        }
                    }
                    
                }        
            }
            
        }
        
       allModels.push_back(*itsCurrentModel); 
       

    }
    
    int  PDBFile::getRecordType(std::string tRecord)
    {   
        for (int i=1; i<(int)PDBRecordHeader.size(); i++)
        {
            Size tCmp = tRecord.compare(0,6,PDBRecordHeader[i]);
            
            if (not tCmp)
            {
                return i;
            }
        }
        
        return 0;
    }
    // Title section
    bool PDBFile::extractTitleInfo(std::string tRecord)
    {
        return true;
    }
    
    // Primary Structure Section
    
    
    bool PDBFile::extractSEQRESInfo(std::string tRecord)
    {
        if(tRecord.size()) 
        {
            
            if( itsCurrentChain->getID() =="")
            {
                itsCurrentChain->setID(tRecord.substr(11,1));
            }
            else if(itsCurrentChain->getID() != tRecord.substr(11,1))
            {
                delete itsCurrentChain;
                itsCurrentChain = NULL;
                itsCurrentChain = new Chain();
                itsCurrentChain->setID(tRecord.substr(11,1));
            }
            
            if (!itsCurrentChain->getNumOfRes())
            {
                itsCurrentChain->setNumOfRes(StrToInt(tRecord.substr(13,4)));
            }
            
            std::vector<std::string> tBuf; 
            StrTokenize(tRecord.substr(19), tBuf);
            
            for (std::vector<std::string>::iterator iT = tBuf.begin();
                    iT != tBuf.end(); ++iT )
            {
                Residue tRes = Residue();
                tRes.setName((*iT));
                itsCurrentChain->residues.push_back(tRes);
            }
            
            if (itsCurrentChain->residues.size()==itsCurrentChain->getNumOfRes())
            {
                allChains.push_back(*itsCurrentChain);
            }
            
            return true;
        }
        else
        {
            return false;
        }
    }
    
    bool PDBFile::extractMODRESInfo(std::string tRecord)
    {
        if(tRecord.size())
        {
            itsCurrentModRes =  new ModRes();
            
            itsCurrentModRes->setID(tRecord.substr(7,4));
            itsCurrentModRes->setName(tRecord.substr(12,3));
            itsCurrentModRes->setChainID(tRecord.substr(16,1));
            itsCurrentModRes->setSeqNum(StrToInt(tRecord.substr(18,4)));
            itsCurrentModRes->setInsCode(tRecord.substr(22,1));   
            itsCurrentModRes->setStdName(tRecord.substr(24,3));
            itsCurrentModRes->setComment(tRecord.substr(29));
        
            allModRes.push_back(*itsCurrentModRes);
            
            return true;
        }
        else 
        {
            return false;
        }
        
    }
    
    
    // Secondary Structure Section
    
    bool PDBFile::extractHelixInfo(std::string tRecord)
    {
        if(tRecord.size())
        {
            Helix * tHelix = new Helix();
            tHelix->setSeriNum(StrToInt(tRecord.substr(7,3)));
            tHelix->setID(tRecord.substr(11,3));
            
            if (!tRecord.substr(15,3).empty())
            {
                Residue * tRes1 = new Residue();
                tRes1->setName(tRecord.substr(15,3));
                tRes1->setChainID(tRecord.substr(19,1));
                tRes1->setSeqNum(StrToInt(tRecord.substr(21,4)));
                tHelix->setInitResidue(*tRes1);
                       
            }
            
            if (!tRecord.substr(27,3).empty())
            {
                Residue * tRes2 = new Residue();
                tRes2->setName(tRecord.substr(27,3));
                tRes2->setChainID(tRecord.substr(31,1));
                tRes2->setSeqNum(StrToInt(tRecord.substr(33,4)));
                tHelix->setEndResidue(*tRes2);        
            }
            
            tHelix->setHelixClass(StrToInt(tRecord.substr(38,2)));
            tHelix->setComment(tRecord.substr(40,30));
            tHelix->setLength(StrToInt(tRecord.substr(71,5)));
            
            allHelices.push_back(*tHelix);
            
            return true;
            
        }
        else 
        {
            return false;
        }
        
    }
    
    bool PDBFile::extractSheetInfo(std::string tRecord)
    {
        if(tRecord.size())
        {
            if (StrToInt(tRecord.substr(7,3)) == 1)
            {   
                itsCurrentSheet = new Sheet();
            }
            
            if ( itsCurrentSheet->getID() =="" )
            {
                itsCurrentSheet->setID(tRecord.substr(11,3));
            }
            if ( itsCurrentSheet->getNumOfStrands() ==0 )
            {
                itsCurrentSheet->setNumOfStrands(StrToInt(tRecord.substr(14,2)));
            }
            // sense of the current strand
            itsCurrentSheet->senses.push_back(StrToInt(tRecord.substr(38,2)));
         
            Strand * tStrand = new Strand();
           
            tStrand->setSeriNum(StrToInt(tRecord.substr(7,3)));
            
            // starting residue in the strand
            Residue * tInitRes=new Residue();
            tInitRes->setName(tRecord.substr(17,3));
            tInitRes->setChainID(tRecord.substr(21,1));
            tInitRes->setSeqNum(StrToInt(tRecord.substr(22,4)));
           
            tStrand->setTerResidue(*tInitRes, 1);
            
            // ending residue in the strand
            Residue * tEndRes = new Residue();
            tEndRes->setName(tRecord.substr(28,3));
            tEndRes->setChainID(tRecord.substr(32,1));
            tEndRes->setSeqNum(StrToInt(tRecord.substr(33,4)));
            
            tStrand->setTerResidue(*tEndRes, -1);
            
            // delete tInitRes and tEndRes
            delete tInitRes;
            tInitRes = NullPoint;
            delete tEndRes;
            tEndRes  = NullPoint;
            
            if (tRecord.size() > 40)
            {
                // Registration atoms
                Atom * tCurAtom = new Atom();
                tCurAtom->setName(tRecord.substr(41,4));
                tCurAtom->setResName(tRecord.substr(45,3));
                tCurAtom->setChainID(tRecord.substr(49,1));
                tCurAtom->setSegNum(StrToInt(tRecord.substr(50,4)));
               
                tStrand->setRegistAtom(*tCurAtom, 1);
               
                Atom * tPrevAtom = new Atom();
                tPrevAtom->setName(tRecord.substr(56,4));
                tPrevAtom->setResName(tRecord.substr(60,3));
                tPrevAtom->setChainID(tRecord.substr(64,1));
                tPrevAtom->setSegNum(StrToInt(tRecord.substr(65,4)));
               
                tStrand->setRegistAtom(*tPrevAtom, -1);
                
                // delete two temp atoms
                delete tCurAtom;
                tCurAtom  = NullPoint;
                delete tPrevAtom;
                tPrevAtom = NullPoint;
                
            }
           
            itsCurrentSheet->allStrands.push_back(*tStrand);

            if( itsCurrentSheet->getNumOfStrands() == tStrand->getSeriNum() ) 
            {
                // the current sheet finishes 
                allSheets.push_back(*itsCurrentSheet);
            }
            
            return true;
                              
        }
        else 
        {
            return false;
        }
    }
    
    // Connectivity Annotation Section
    
    bool PDBFile::extractSSBondInfo(std::string tRecord)
    {
        if(tRecord.size())
        {
            SSBond * tSSBond = new SSBond();
           
            tSSBond->setName(tRecord.substr(0,6));
            tSSBond->setSeriNum(StrToInt(tRecord.substr(7,3)));
            
            Residue * tResidue1 = new Residue();
            tResidue1->setName(tRecord.substr(11,3));
            tResidue1->setChainID(tRecord.substr(15,1));
            int tI = StrToInt(tRecord.substr(17,4));
            tResidue1->setSeqNum(tI);
          
            tSSBond->resSeqNums[0] =tI;
            
            Residue * tResidue2 = new Residue();
            tResidue2->setName(tRecord.substr(26,3));
            tResidue2->setChainID(tRecord.substr(29,1));
            tI = StrToInt(tRecord.substr(31,4));
            tResidue2->setSeqNum(tI);           
            
            tSSBond->resSeqNums[1] = tI;
            
            tSSBond->residues.push_back(*tResidue1);
            tSSBond->residues.push_back(*tResidue2);
            
            tSSBond->resSym[0]  = StrToInt(tRecord.substr(59,6));
            tSSBond->resSym[1]  = StrToInt(tRecord.substr(66,6));  
            
            tSSBond->setLength(StrToReal(tRecord.substr(73,5)));
            
            allSSBonds.push_back(*tSSBond);
            
        }
        
        return true;
    }
    
    bool PDBFile::extractLinkInfo(std::string tRecord)
    {
        if(tRecord.size())
        {
            Link * tLink = new Link();
            
            tLink->setName(tRecord.substr(0,6));
            
            Atom * tAtom1 = new Atom();
            tAtom1->setName(TrimSpaces(tRecord.substr(12,4)));
            tAtom1->setAltLoc(tRecord.substr(16,1));
            tAtom1->setResName(tRecord.substr(17,3));
            tAtom1->setChainID(tRecord.substr(21,1));
            tAtom1->setSegNum(StrToInt(tRecord.substr(22,4)));
            tAtom1->setInsCode(tRecord.substr(26,1));
            
            tLink->atoms.push_back(*tAtom1);
            delete tAtom1;
            tAtom1 = NULL;
            
            Atom * tAtom2 = new Atom();
            tAtom2->setName(TrimSpaces(tRecord.substr(42,4)));
            tAtom2->setAltLoc(tRecord.substr(46,1));
            tAtom2->setResName(tRecord.substr(47,3));
            tAtom2->setChainID(tRecord.substr(51,1));
            tAtom2->setSegNum(StrToInt(tRecord.substr(52,4)));
            tAtom2->setInsCode(tRecord.substr(56,1));
            
            tLink->atoms.push_back(*tAtom2);
            delete tAtom2;
            tAtom2 = NULL;
            
            tLink->symOp[0] = StrToInt(tRecord.substr(59,6));
            
            tLink->symOp[0] = StrToInt(tRecord.substr(66,6));
            
            tLink->setLength(StrToReal(tRecord.substr(73,5)));            
            
            allLinks.push_back(*tLink);
            
            delete tLink;
            tLink = NULL;
            
        }
        
        return true;
    }
    
    // Crystallographic and Coordinate Transformation Section
    bool PDBFile::extractCryst1Info(std::string tRecord)
    {
        if(tRecord.size())
        {
            /*
            itsCrystInfo=new CrystInfo(StrToReal(tRecord.substr(6,9)),
                    StrToReal(tRecord.substr(15,9)),
                    StrToReal(tRecord.substr(24,9)),
                    StrToReal(tRecord.substr(33,7)),
                    StrToReal(tRecord.substr(40,7)),
                    StrToReal(tRecord.substr(47,7)),
                    tRecord.substr(55,11),
                    StrToInt(tRecord.substr(66,4)));
           */
            return true;
        }
        
        return false;        
    }
    

    bool PDBFile::extractOrigXNInfo(std::string tRecord)
    {
        if (tRecord.size())
        {
            itsCrystInfo->ORIGXn.push_back(StrToReal(tRecord.substr(10,10)));
            itsCrystInfo->ORIGXn.push_back(StrToReal(tRecord.substr(20,10)));
            itsCrystInfo->ORIGXn.push_back(StrToReal(tRecord.substr(30,10)));
            itsCrystInfo->ORIGXn.push_back(StrToReal(tRecord.substr(45,10)));
         
            return true;
        }
         
        return false;
        
    }
 
    bool PDBFile::extractMatrixNInfo(std::string tRecord)
    {
        if (tRecord.size())
        {
            itsCrystInfo->MTRIXn.push_back(StrToReal(tRecord.substr(10,10)));
            itsCrystInfo->MTRIXn.push_back(StrToReal(tRecord.substr(20,10)));
            itsCrystInfo->MTRIXn.push_back(StrToReal(tRecord.substr(30,10)));
            itsCrystInfo->MTRIXn.push_back(StrToReal(tRecord.substr(45,10)));
         
            return true;
        }
         
        return false;
        
    }
    
    bool PDBFile::extractScaleNInfo(std::string tRecord)
    {
        if (tRecord.size())
        {
            itsCrystInfo->SCALEn.push_back(StrToReal(tRecord.substr(10,10)));
            itsCrystInfo->SCALEn.push_back(StrToReal(tRecord.substr(20,10)));
            itsCrystInfo->SCALEn.push_back(StrToReal(tRecord.substr(30,10)));
            itsCrystInfo->SCALEn.push_back(StrToReal(tRecord.substr(45,10)));
         
            return true;
        }
         
        return false;
        
    }    
    
    
    // Coordinate Section
    void PDBFile::initOneModel(std::string tRecord)
    {
        if (itsCurrentModel == NULL)
        {
            itsCurrentModel = new Model();
            itsCurrentModel->setSeriNum(1);
        }
    }
    
    bool PDBFile::extractAtomInfo(std::string tRecord)
    {
        
        /*
        if (itsCurrentModel == NULL)
        {
            itsCurrentModel = new Model();
            itsCurrentModel->setSeriNum(allModels.size()+1);
        }
        if (itsCurrentChain == NULL)
        {
            itsCurrentChain = new Chain();
            itsCurrentChain->setModSeriNum(itsCurrentModel->getSeriNum());
        }
        
        if (itsCurrentResidue == NULL)
        {
            itsCurrentResidue = new Residue();
        }
        */
        if(!itsCurrentAltLoc)
        {
        Atom  tAtom;
        
        tAtom.setModSeriNum(itsCurrentModel->getSeriNum());
        //std::cout << "Atom's model numb : " << tAtom.getModSeriNum() << std::endl;
        tAtom.setSeriNum(StrToInt(tRecord.substr(6,5)));
        //std::cout << "Atom's serial numb : " << tAtom.getSeriNum() << std::endl;
        tAtom.setName(TrimSpaces(tRecord.substr(12,4)));
        //std::cout << "Atom's Name: " << tAtom.getName() << std::endl;
        tAtom.setAltLoc(tRecord.substr(16,1));
        if(TrimSpaces(tRecord.substr(16,1)).length() !=0)
        {
            itsCurrentAltLoc = true;
        }
        //std::cout << "Atom's AltLoc: " << tAtom.getAltLoc() << std::endl;
        Name tResName = tRecord.substr(17,3);
        tAtom.setResName(TrimSpaces(tResName));
        //std::cout << "Atom's Residue Name: " << tAtom.getResName() << std::endl;
        Name tCName=tRecord.substr(21,1);
        
        tAtom.setChainID(TrimSpaces(tCName));
        //std::cout << "Atom'Chain ID: " << tAtom.getChainID() << std::endl;
        tAtom.setSeqNum(StrToInt(tRecord.substr(22,4)));
        //std::cout << "Atom's residue seq numb: " << tAtom.getSeqNum() << std::endl;
        if(TrimSpaces(tRecord.substr(26,1)).length() !=0)
        {
            tAtom.setInsCode(TrimSpaces(tRecord.substr(26,1)));
        }
        
        //std::cout << "Atom's InsCode: " << tAtom.getInsCode() << std::endl;
        tAtom.coords.push_back(StrToReal(tRecord.substr(30,8)));
        tAtom.coords.push_back(StrToReal(tRecord.substr(38,8))); 
        tAtom.coords.push_back(StrToReal(tRecord.substr(46,8)));
        // std::cout << "Atom's coord z: " << tAtom.coords[2] << std::endl;
        
        tAtom.setOccup(StrToReal(tRecord.substr(54,6)));
        tAtom.setTempFact(StrToReal(tRecord.substr(60,6)));
        
        tAtom.setElementType(TrimSpaces(tRecord.substr(76,2)));
        //std::cout << "Atom's Element type : " 
        //          << tAtom.getElementType() << std::endl;
        tAtom.setPartialCharge(StrToReal(tRecord.substr(78,2)));
        
        
        allAtomList.push_back(tAtom);
     
        //int ii = allAtomList.size() -1;
        //std::cout<< "atom name            : " << allAtomList[ii].getName() 
        //        << std::endl;
        //std::cout<< "atom seriNum        : " << allAtomList[ii].getSeriNum() 
        //        << std::endl;
       
        //std::cout<< "atom Residue Name   : " << allAtomList[ii].getResName()
        //        << std::endl;
        //std::cout<< "atom ChainID        : " << allAtomList[ii].getChainID() 
        //        << std::endl;
        //std::cout << "Atoms resSeq numb  : " << allAtomList[ii].getSeqNum() 
        //        << std::endl;
        //std::cout << "Atom  coords Z     : " << allAtomList[ii].coords[2] 
        //          << std::endl;
        //std::cout << "Atom Element type : " 
        //          << allAtomList[ii].getElementType() << std::endl;
    
        
        // add it to the Het-atom list if it is
     
        if (!tRecord.substr(0,6).compare("HETATM"))
        {
            allHetAtmList.push_back(tAtom);
        }

     
     
        if (itsCurrentResidue->getName() == NullString)
        {
            itsCurrentResidue->setName(tAtom.getResName());
        }
       
        if (itsCurrentResidue->getSeqNum() == ZeroInt)
        {
           itsCurrentResidue->setSeqNum(tAtom.getSeqNum());
        } 
        
        if (itsCurrentResidue->getChainID() == NullString)
        {
           itsCurrentResidue->setChainID(tAtom.getChainID());
        }
        
        if (itsCurrentResidue->getInsCode() == NullString)
        {
           itsCurrentResidue->setInsCode(tAtom.getInsCode());
        } 
        
        // add the atom to a residue
        //std::cout << tAtom.getResName() << std::endl;
        //std::cout << itsCurrentResidue->getName() << std::endl;
        //std::cout << tAtom.getSeqNum() << std::endl;
        //std::cout << itsCurrentResidue->getName() << std::endl;
          
        if (tAtom.getSeqNum() == itsCurrentResidue->getSeqNum())
        {
           itsCurrentResidue->atoms.push_back(tAtom);
           // std::cout <<  itsCurrentResidue->atoms.size() << std::endl;
        }
        else 
        {
            // add the current residue to a chain

            itsCurrentChain->addOneResidue(*itsCurrentResidue);
           
            int iR = itsCurrentChain->residues.size()-1;
            itsCurrentChain->residues[iR].atoms.clear();
            for (int i = 0; i < (int)itsCurrentResidue->atoms.size(); i++)
            {
                itsCurrentChain->residues[iR].atoms.push_back(itsCurrentResidue->atoms[i]);
            }
           
            if (itsCurrentChain->getID() == "")
            {
                itsCurrentChain->setID(tAtom.getChainID());
            }
            
            itsCurrentResidue = new Residue();
            itsCurrentResidue->setName(tAtom.getResName());
            itsCurrentResidue->setSeqNum(tAtom.getSeqNum());
            itsCurrentResidue->setChainID(tAtom.getChainID());
            itsCurrentResidue->setModelSeriNum(tAtom.getModSeriNum());
            itsCurrentResidue->atoms.push_back(tAtom);
        }
        
        }
        else
        {
            if (tRecord.substr(16,1) !="")
            {
                allAtomList.back().altLoc.push_back(StrToReal(tRecord.substr(30,8)));
                allAtomList.back().altLoc.push_back(StrToReal(tRecord.substr(38,8)));
                allAtomList.back().altLoc.push_back(StrToReal(tRecord.substr(46,8)));
                
                itsCurrentResidue->atoms.back().altLoc.push_back(StrToReal(tRecord.substr(30,8)));
                itsCurrentResidue->atoms.back().altLoc.push_back(StrToReal(tRecord.substr(38,8)));
                itsCurrentResidue->atoms.back().altLoc.push_back(StrToReal(tRecord.substr(46,8)));
                itsCurrentAltLoc = false;
            }
        }
        
        return true;
    }
    
    bool PDBFile::extractAnisouInfo(std::string tRecord)
    {
        itsCurrentAtom->Uxx.push_back(StrToInt(tRecord.substr(28,7)));
        itsCurrentAtom->Uxx.push_back(StrToInt(tRecord.substr(35,7)));
        itsCurrentAtom->Uxx.push_back(StrToInt(tRecord.substr(42,7)));
        itsCurrentAtom->Uxx.push_back(StrToInt(tRecord.substr(49,7)));
        itsCurrentAtom->Uxx.push_back(StrToInt(tRecord.substr(56,7)));
        itsCurrentAtom->Uxx.push_back(StrToInt(tRecord.substr(63,7)));
        return true;
    }
    
    // tRecord is not used yet
    bool PDBFile::extractTerInfo(std::string tRecord)
    {     
        itsCurrentChain->addOneResidue(*itsCurrentResidue);
        itsCurrentModel->addOneChain(*itsCurrentChain);
        
        delete itsCurrentResidue;
        itsCurrentResidue = NULL;
        itsCurrentResidue = new Residue();
        
        delete itsCurrentChain;
        itsCurrentChain = NULL;
        itsCurrentChain = new Chain();
        return true;
    }
    
    // tRecord is not used yet
    bool PDBFile::extractEndMdlInfo(std::string tRecord)
    {
        allModels.push_back(*itsCurrentModel);
        return true;
    }
    
    bool PDBFile::extractConnectInfo(std::string tRecord)
    {
        int tSeriNum = StrToInt(tRecord.substr(0,6));
        
        for (std::vector<Atom>::iterator ia=allAtomList.begin();
                ia !=allAtomList.end(); ia++)
            if(tSeriNum==ia->getSeriNum())
            {
                
                ia->connectedAtoms.push_back(tSeriNum);
                return true;
            }
        
        return false;     
    }
    
    // get and set methods
    ID PDBFile::getPDBID() const
    {
        return allHeaderInfo.HEADER->PDBId;
    }
    
    Date PDBFile::getDepositedDate() const
    {
        return allHeaderInfo.HEADER->depositDate;
    }
        
       
    DictPDBFile::DictPDBFile() 
    {
    }
    
    DictPDBFile::DictPDBFile(Name tFname, std::ios_base::openmode tOpenMode)
    {
        if (tOpenMode == std::ios::in)
        {
            inFile.open(tFname.c_str(), tOpenMode);
            setupSystem();
        }
    }
    
    void DictPDBFile::setupSystem()
    {
        if (inFile.is_open() )
        { 
            std::string tRecord="";
            while(!inFile.eof())
            {
                std::getline(inFile, tRecord);
                std::vector<std::string> tStrs;
                StrTokenize(tRecord, tStrs);
                if ((int)tStrs.size() >=12 )
                {
                    if (tStrs[0].find("HETATM") !=std::string::npos 
                         || tStrs[0].find("ATOM") !=std::string::npos)
                    {
                        ID aID(TrimSpaces(tRecord.substr(12,4)));
                        allHetAtmList[aID].push_back(StrToReal(tRecord.substr(30,8)));
                        allHetAtmList[aID].push_back(StrToReal(tRecord.substr(38,8)));
                        allHetAtmList[aID].push_back(StrToReal(tRecord.substr(46,8)));
                    }
                }
            }
            inFile.close();
        }
    }
    
    void extern outPDB(FileName tFName, 
                       ID tMonoRootName,
                       std::vector<LIBMOL::AtomDict>& tAtoms)
    {
        // This is a temporary one, the method should be defined outside 
        // this class.
        std::string tName(tFName);    
       
        std::vector<std::string> parts;
        StrTokenize(tName, parts, '.');
        std::string outPDBName = parts[0] + ".pdb";
        
        std::ofstream outPDB(outPDBName.c_str());
        
        if(outPDB.is_open())
        {
            // Header section
            
            srand((unsigned)std::time( NULL ));
            outPDB.width(10);
            outPDB << std::left << "HEADER";
            outPDB.width(30);
            outPDB << std::left << " MONOMER tests " << tMonoRootName;
            outPDB.width(10);
            outPDB << std::left << " Date:";
            outPDB.width(6);
            outPDB << std::left << tMonoRootName;
            outPDB.width(14);
            outPDB << std::left <<"" <<std::endl;
                
            // CRYST1 section 
            outPDB.width(80);
            outPDB << std::left 
                   <<"CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1" 
                   <<std::endl;
            
            // ATOM sections
            for (std::vector<AtomDict>::iterator iA=tAtoms.begin();
                    iA !=tAtoms.end(); iA++)
            {
                
                std::string tID(iA->id);
                if (tID.find("\"") !=tID.npos)
                {
                    char tQ='\"';
                    cleanChar(tID, tQ);
                }
                
                //double r1 =  (double) rand()/RAND_MAX;
                //double r2 =  (double) rand()/RAND_MAX;
                //double r3 =  (double) rand()/RAND_MAX;
                
                outPDB.width(6);
                outPDB <<std::left<< "HETATM";
                outPDB.width(5);
                outPDB <<std::right << iA->seriNum+1;
                outPDB.width(1);
                outPDB << " ";
                
                if ((int)tID.size() <=4)
                {
                    outPDB.width(4);
                    outPDB << std::left << tID;
                }
                else
                {
                    outPDB.width(4);
                    outPDB << std::left << tID.substr(0,4);
                }
                
                outPDB.width(1); // altLoc
                outPDB << " ";
                outPDB.width(3); // resName
                outPDB << tMonoRootName.substr(0,3);
                outPDB.width(1); // empty
                outPDB << " ";
                outPDB.width(1);  // chainID 
                outPDB << std::right << "A";
                outPDB.width(4);  // resSeq
                outPDB << std::right << "1";
                outPDB.width(1);  // iCode
                outPDB << " "; 
                outPDB.width(3); // empty
                outPDB << "  "; 
                outPDB.width(8);
                outPDB << std::right << std::setprecision(3) 
                        <<std::fixed << iA->coords[0];
                outPDB.width(8);
                outPDB << std::right << std::setprecision(3) 
                        <<std::fixed << iA->coords[1];
                outPDB.width(8);
                outPDB << std::right << std::setprecision(3) 
                        <<std::fixed << iA->coords[2];
                
                outPDB.width(6);
                outPDB << std::right << "1.00";
                outPDB.width(6);
                outPDB << std::right << std::setprecision(2) << std::fixed << "20.00";
                outPDB.width(12);
                outPDB << std::right << iA->chemType;
                outPDB.width(2);
                outPDB << std::right << "" << std::endl;
            }
        }
        
        outPDB.close();
    }
}
