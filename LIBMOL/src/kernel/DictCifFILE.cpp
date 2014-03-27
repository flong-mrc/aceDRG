/*
 * 
 * File:   CIFFile.cpp
 * Author: flong
 *
 * last updated on  Jan 24, 2012, 05:36 PM
*/

#include "DictCifFile.h"

namespace LIBMOL
{
    
    GenCifFile::GenCifFile() : curBlockLine(ZeroInt),
            hasCoords(false),
            hasH(false),
            itsCurAtomSeriNum(ZeroInt),
            itsCurAtom(NullPoint),
            itsCurBlock(""),
            itsCurCryst(NullPoint)
    {
    }
    
    GenCifFile::GenCifFile(Name tFname, 
                           std::ios_base::openmode tOpenMode) :
                           curBlockLine(ZeroInt),
                           hasCoords(false),
                           hasH(false),
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
                setupSystem();
            }
            else
            {
                std::cout << tFname << " Can not be opened for reading " << std::endl;
                exit(1);
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
                           hasCoords(false),
                           hasH(false),
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
            
            std::vector<std::string>     tBlocLines;
            
            while(!inFile.eof())
            {   
                std::getline(inFile, tRecord);
                tRecord = TrimSpaces(tRecord);
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
            
            
            
            /*
            std::cout << "There are " << (int)allHydroAtoms.size() 
                    << " Hydrogen atoms in the system. They are " << std::endl;
            
            for (std::vector<int>::iterator iA=allHydroAtoms.begin();
                   iA !=allHydroAtoms.end(); iA++)
            {
                std::cout << "Atom " << allAtoms[*iA].id << std::endl;
            }
          
            std::cout << "There are " << (int)allBonds.size() 
                    << " bonds in the system. They are : " << std::endl;
            for (std::vector<BondDict>::iterator iB = allBonds.begin();
                    iB !=allBonds.end(); iB++)
            {
                std::cout << "Bond between atom " << iB->atoms[0]  
                          << " and atom " << iB->atoms[1] << std::endl ; 
            }
            
            std::cout << "There are " << (int)allChirals.size() 
                      << " chirals in the system. They are: " <<std::endl;
            
            for (std::vector<ChiralDict>::iterator iCh = allChirals.begin();
                        iCh != allChirals.end(); iCh++)
            {
                std::cout << iCh->id  << "\t";
                for (std::vector<int>::iterator iAt=iCh->atoms.begin();
                            iAt != iCh->atoms.end(); iAt++)
                {
                    std::cout << allAtoms[*iAt].id << "\t";
                }
                std::cout << iCh->sign << std::endl;
            }
                
            setHydroAtomConnect();
            
            // addMissHydroAtoms();
            
            setAtomsBondingAndChiralCenter();
            
            // setAllAngles();
            
            setAtomsCChemType();
            
            setAtomsMetalType();
            
            setAtomsVDWRadius();
            
            setAtomsPartialCharges();
            
            ringDetecting();
            
            setAtomsCCP4Type();
        
           
            setAllAngles();
            ringDetecting();
            
            for (std::map<ID, std::vector<RingDict> >::iterator iMR=allRings.begin();
                    iMR != allRings.end(); iMR++)
            {
                for (std::vector<RingDict>::iterator iR=iMR->second.begin();
                        iR != iMR->second.end(); iR++)
                {
                    allRingsV.push_back(*iR);
                }
            }
            detectPlaneGroups();
            
            //setAllTorsions();
            setAllTorsions2();
            */
            
            //buildAtomTree tTool;
            
            //tTool.buildTree(allAtoms, allBonds, allAngles, allTorsions, 
            //                allRingsV, allPlanes, allChirals);
            
           
            
            
            // std::cout << "finish reading the input cif file: " 
            //          << std::endl
            /*
            if (hasCoords)
            {
                std::cout << "The system has coords " << std::endl;
            }
            else
            {
                std::cout << "The system has no coords " << std::endl;
            }
            for (std::vector<AtomDict>::iterator iA = allAtoms.begin();
                    iA != allAtoms.end(); iA++)
            {
                std::cout << "\nAtom " << iA->seriNum << " : " << std::endl
                        << "Its ID : " << iA->id << std::endl
                        << "Its Chemical Type : " << iA->chemType << std::endl
                        << "Its COD chemical type " << iA->cChemType << std::endl
                        << "Its bonding index : "   << iA->bondingIdx << std::endl
                        << "Its residue Name: " << iA->resName<< std::endl;
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
             */
        }
        
       
        
    } 
    
    
    void GenCifFile::initAllCifKeys()
    {
        std::string clibMonDir(std::getenv("CLIBD_MON"));
        std::string fName(clibMonDir);
        fName.append("/list/cif_tag.list");
        //std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
        //std::string fName(clibMonDir);
        //fName.append("/lib/cif_tag.list");
        
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
                    std::cout << "Wrong Line format for the lattice in line " << *iBl << std::endl;
                    exit(1);
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
                    exit(1);
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
                    exit(1);
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
                    exit(1);
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
                    exit(1);
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
                {
                    std::cout << "No cell parameter c in line " << *iBl << std::endl;
                    exit(1);
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
                        StrTokenize(tBuf[1], tBuf2); 
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
                    exit(1);
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
                        StrTokenize(tBuf[1], tBuf2); 
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
                    exit(1);
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
                        StrTokenize(tBuf[1], tBuf2); 
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
                    exit(1);
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
                    std::cout << "No cell parameter c in line " << *iBl << std::endl;
                    exit(1);
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
                    //std::cout << curBlockLine << std::endl;
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
                std::cout << *iBl << std::endl;
                
                itsCurAtom = new AtomDict();
                bool getCoords = false;
                itsCurAtom->seriNum = itsCurAtomSeriNum;
                itsCurAtomSeriNum ++;
                if (hasProps["atom"].find("id") != hasProps["atom"].end())
                {
                    itsCurAtom->existProps["id"] =  hasProps["atom"]["resName"];
                    cleanSymbol(tF[hasProps["atom"]["id"]], "\"");
                    itsCurAtom->id = TrimSpaces(tF[hasProps["atom"]["id"]]);
                    std::cout << "Its ID: " << itsCurAtom->id << std::endl;
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
                    std::cout << "Its chemType : " << itsCurAtom->chemType << std::endl;
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
                    //std::cout << "Its U_iso : " << itsCurAtom->isoB << std::endl;
                }
                if (hasProps["atom"].find("parCharge") != hasProps["atom"].end())
                {
                    itsCurAtom->existProps["parChange"] = hasProps["atom"]["parCharge"];
                    itsCurAtom->parCharge = StrToReal(tF[itsCurAtom->existProps["parChange"]]);
                    //std::cout << "Its partialCharge :" 
                    //        << itsCurAtom->parCharge
                    //        << std::endl;
                }
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
                   
                    std::cout << "Its (fractional) coord x : " << itsCurAtom->fracCoords[0] << std::endl;
                    
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
                            std::cout << "the string for coordinate y is " << tSY << std::endl;
                        }
                    }
                    std::cout << "Its (fractional) coord y : " << itsCurAtom->fracCoords[1] << std::endl;
                    
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
                            std::cout << "the string for coordinate z is " << tSZ << std::endl;
                        }
                    }
                    std::cout << "Its (fractional) coord z : " << itsCurAtom->fracCoords[2] << std::endl;
                    
                    getCoords = true;
                    
                }
               
                if (itsCurAtom->existProps.find("fract_x") !=itsCurAtom->existProps.end()
                    && getCoords)
                {
                    // TranslateIntoUnitCell(itsCurAtom->fracCoords);
                    FractToOrtho(itsCurAtom->fracCoords, itsCurAtom->coords, itsCurCryst->itsCell->a, 
                            itsCurCryst->itsCell->b, itsCurCryst->itsCell->c, itsCurCryst->itsCell->alpha,
                            itsCurCryst->itsCell->beta, itsCurCryst->itsCell->gamma);
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
    
   
    /*
    void GenCifFile::getCifAtomInfo(std::vector<std::vector<std::string> >::iterator iBs)
    {
        
    }
     */
    
    // ################################## another class for cif files of CCP4 
    // dictionary 
    
    
    
    DictCifFile::DictCifFile() : curBlockLine(ZeroInt),
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
                        //std::cout << "get compound info" << std::endl;
                        getChemInfo(tRecord);
                    }
                    else if (lBloc["dataDesc"])
                    {
                        getDataDescription(tBuf);
                    }
                    else if(lBloc["atom"])
                    {   
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
            
            // Set the bonding properties for all atoms based
            // on their connections
            // If no bonds defined in the cif file. Stop the program at the moment.
            if ((int)allBonds.size()==0)
            {
                std::cout << "There is no bond defined in the cif file. Program stopped"
                        << std::endl;
                exit(1);
            }
            
            
                
            setHydroAtomConnect();
            // addMissHydroAtoms();
            
            setAtomsBondingAndChiralCenter();
            
            // setAllAngles();
            
            setAtomsCChemType();
            
            setAtomsMetalType();
            
            setAtomsVDWRadius();
            
            setAtomsPartialCharges();
            
            ringDetecting();
            
            if (!hasCCP4Type)
            {
                setAtomsCCP4Type();
            }
            
            /*
            std::map<unsigned, ID> tCMap;
            // check 
            for (unsigned i=0; i < allAtoms.size(); i++)
            {
                tCMap[i] = allAtoms[i].ccp4Type;
            }
            setAtomsCCP4Type();
            
            std::cout << "CCP4 atom types : " << std::endl;
            
            for (unsigned i=0; i < allAtoms.size(); i++)
            {
                if (allAtoms[i].ccp4Type.compare(tCMap[i]) !=0)
                {
                    std::cout << "For atom " << i << " id : " << allAtoms[i].id
                            << std::endl << "Its CCP4 type: " << allAtoms[i].ccp4Type
                            << " and " << tCMap[i] << std::endl;
                }
            }
                
            */
            //}
          
            /*
            setAllAngles();
            ringDetecting();
            
            for (std::map<ID, std::vector<RingDict> >::iterator iMR=allRings.begin();
                    iMR != allRings.end(); iMR++)
            {
                for (std::vector<RingDict>::iterator iR=iMR->second.begin();
                        iR != iMR->second.end(); iR++)
                {
                    allRingsV.push_back(*iR);
                }
            }
            detectPlaneGroups();
            
            //setAllTorsions();
            setAllTorsions2();
            */
            //buildAtomTree tTool;
            
            //tTool.buildTree(allAtoms, allBonds, allAngles, allTorsions, 
            //                allRingsV, allPlanes, allChirals);
            
           
            
            // std::cout << "finish reading the input cif file: " 
            //          << std::endl
            
            if (hasCoords)
            {
                std::cout << "The system has coords " << std::endl;
            }
            else
            {
                std::cout << "The system has no coords " << std::endl;
            }
            for (std::vector<AtomDict>::iterator iA = allAtoms.begin();
                    iA != allAtoms.end(); iA++)
            {
                std::cout << "\nAtom " << iA->seriNum << " : " << std::endl
                        << "Its ID : " << iA->id << std::endl
                        << "Its Chemical Type : " << iA->chemType << std::endl
                        << "Its CCP4 chemical type " << iA->ccp4Type << std::endl
                        << "Its bonding index : "   << iA->bondingIdx << std::endl
                        << "Its residue Name: " << iA->resName<< std::endl;
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
                //std::cout << "its order : " << iBo->order << std::endl;
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
                else if (TrimSpaces(tF[0]).find("_chem_comp.") 
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
                    std::cout<< "Setup flag atom " <<std::endl;
                    curBlockLine = 0;
                    setFlags(tL, "atom");
                }
                else if ( ! tL["bond"] && 
                         TrimSpaces(tF[0]).find("_chem_comp_bond")
                         !=std::string::npos)
                {
                    std::cout << "setup bond flag " << std::endl;
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
           
           if(tF[0].find("_name") !=std::string::npos)
           {
               dictCifHead.libName = TrimSpaces(tF[1]);
           }
           else if (tF[0].find("_version") !=std::string::npos)
           {
               dictCifHead.libVersion = TrimSpaces(tF[1]);
           }
           else if (tF[0].find("_update") !=std::string::npos)
           {
               dictCifHead.libUpdate = TrimSpaces(tF[1]);
           }
        }
        
        //std::cout << "dictCifHead.libName: " <<  dictCifHead.libName << std::endl;
        //std::cout << "dictCifHead.libUpdate: " <<  dictCifHead.libUpdate << std::endl;
        //std::cout << "dictCifHead.libVersion: " <<  dictCifHead.libVersion << std::endl;
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
                tF4.push_back(tF[1]);
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
                        // std::cout << "Compound name " << propComp.name  << std::endl;
                    }
                    if (hasProps["compoundInfo"].find("group") != hasProps["compoundInfo"].end())
                    {
                        propComp.group = TrimSpaces(tF4[hasProps["compoundInfo"]["group"]]);
                        //std::cout << "Compound group " << propComp.group  << std::endl;
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
                            // std::cout << "number of H atoms " << propComp.numH  << std::endl;
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
                     tID) != allCifKeys["dataDesc"].end()); 
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
                    //std::cout << curBlockLine << std::endl;
                    //curBlockLine++;
                }
                else if(tF1[1].find("charge") !=std::string::npos)
                {
                    hasProps["atom"].insert(std::pair<std::string, int>("parCharge",curBlockLine));
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
        
        //if ((int)tF.size() == (int)hasProps["atom"].size() 
        //        && (int)tF.size() >2 && tF[0].find("#") ==std::string::npos)
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
                // std::cout << "Its ID: " << itsCurAtom->id << std::endl;
            }
            if (hasProps["atom"].find("chemType") != hasProps["atom"].end())
            {
                itsCurAtom->existProps["chemType"] =  hasProps["atom"]["chemType"];
                itsCurAtom->chemType = TrimSpaces(tF[hasProps["atom"]["chemType"]]);
                
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
                //std::cout << "Its partialCharge :" 
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
                        std::cout << "Its coord y : " << itsCurAtom->coords[1] << std::endl;
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
                itsCurBond->order = tF[hasProps["bond"]["order"]];
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
                    allAtoms[iPos2].connAtoms.push_back(iPos1);
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
    
    void DictCifFile::setAtomsBondingAndChiralCenter()
    {
        // First round
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {
            int t_len =0;
            for (std::vector<int>::iterator iConn=iAt->connAtoms.begin();
                    iConn !=iAt->connAtoms.end(); iConn++)
            {
                if(!allAtoms[*iConn].isMetal)
                {
                    t_len++;
                }
            }

            if (iAt->chemType.compare("C")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    if (iAt->chiralIdx ==0)
                    {
                        iAt->chiralIdx  = 2;
                    }
                    iAt->bondingIdx = 3;
                }
                else if (t_len ==3)
                {
                    iAt->chiralIdx  = 0;
                    iAt->bondingIdx = 2;
                } 
            }
            else if (iAt->chemType.compare("N")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    if (iAt->chiralIdx ==0)
                    {
                        iAt->chiralIdx  = 2;
                    } 
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3;  
                }
                //else if (t_len==3) // temp
                //{
                //    iAt->chiralIdx  = -1;
                //    iAt->bondingIdx =  2;
                //}
                else if (t_len ==2)
                {
                    iAt->chiralIdx  = 0;
                    iAt->bondingIdx = 2;
                } 
            }
            else if (iAt->chemType.compare("B")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3;
                }
            }
            else if (iAt->chemType.compare("O")==0)
            {
                if ((int)iAt->connAtoms.size()==2)
                {
                    iAt->bondingIdx = 2;
                }
            }
            else if (iAt->chemType.compare("SI")==0 
                    || iAt->chemType.compare("P")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    if (iAt->chiralIdx ==0)
                    {
                        std::vector<ID> atps;
                        for (std::vector<int>::iterator iNA=iAt->connAtoms.begin();
                                iNA !=iAt->connAtoms.end(); iNA++)
                        {
                            if (std::find(atps.begin(), atps.end(), allAtoms[*iNA].chemType)==atps.end())
                            {
                                atps.push_back(allAtoms[*iNA].chemType);
                            }
                        }
                        if ((int)atps.size() >2)
                        {
                            iAt->chiralIdx  = 2;
                        }
                        else
                        {
                            iAt->chiralIdx =0;
                        }
                    }
                   
                    iAt->bondingIdx = 3; 
                }
                else if (t_len==3)
                {
                    if (iAt->chiralIdx ==0)
                    {
                        iAt->chiralIdx  = 2;
                    }
                    
                    iAt->bondingIdx = 2; 
                }
            }
            else if (iAt->chemType.compare("S")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4 || t_len==3)
                {
                    if (iAt->chiralIdx ==0)
                    {
                        iAt->chiralIdx  = 2;
                    }
                    iAt->bondingIdx = 3; 
                }
            }
        }
        
        // more conditions 
        
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {   
            int t_len =0;
            for (std::vector<int>::iterator iConn=iAt->connAtoms.begin();
                    iConn !=iAt->connAtoms.end(); iConn++)
            {
                if(!allAtoms[*iConn].isMetal)
                {
                    t_len++;
                }
            }

            if (iAt->chemType.compare("N")==0 || iAt->chemType.compare("B")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==3)
                {
                    bool l_sp2 = false;
                    for (std::vector<int>::iterator iCA=iAt->connAtoms.begin();
                            iCA != iAt->connAtoms.end(); iCA++)
                    {
                        if(allAtoms[*iCA].bondingIdx == 2)
                        {
                            l_sp2 = true;
                        }
                    }
                    if (l_sp2)
                    {
                        // Now we can say this atom is in sp2 orbits 
                        iAt->chiralIdx  =  0;
                        iAt->bondingIdx =  2;
                    }
                    else
                    {
                        if (iAt->chiralIdx ==0)
                        {
                            iAt->chiralIdx  = 2;
                        }
           
                        iAt->bondingIdx =  3;
                    }
                } 
            }
        }
        
        // Further check if a chiral center is a real one
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA != allAtoms.end(); iA++)
        {
            if (iA->chiralIdx !=0)
            {
                // std::cout << "Atom " << iA->id << std::endl;
                        
                std::vector<ID> chirRAtms;
                for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                        iNB != iA->connAtoms.end(); iNB++)
                {
                    //std::cout << "NB atom " <<allAtoms[*iNB].chemType << std::endl;
                    std::size_t tFind = allAtoms[*iNB].chemType.find("H");
                    if (tFind !=std::string::npos)
                    {
                        chirRAtms.push_back(allAtoms[*iNB].id);
                    }
                }
                // std::cout << "H atom " << (int)chirRAtms.size() << std::endl;
                if ((int)chirRAtms.size() >1 && (int)iA->connAtoms.size() <=4 )
                {
                    iA->chiralIdx = 0;
                }
                //std::cout << "its chiral idx " << iA->chiralIdx << std::endl;
            }
        }
        
        // No need for the third round, those could be defined in 
        // the first round
        // Check
        
        std::cout << "Chiral and plane feather for atoms in the system" 
                  << std::endl;
        
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {
            std::cout << "Atom " << iAt->id << " with " 
                    << (int)iAt->connAtoms.size() << " connected." 
                    << std::endl;
            if (iAt->bondingIdx == 1)
            {
                std::cout << "Atom " << iAt->id << " may be in lines " 
                        << std::endl;
            }
            else if (iAt->bondingIdx == 2)
            {
                std::cout << "Atom " << iAt->id << " may be in planes " 
                        << std::endl;
            }
            else if (iAt->bondingIdx == 3)
            {
                std::cout << "Atom " << iAt->id << " is at a tetrahedra center " 
                        << std::endl;
            }
            else 
            {
                std::cout << "Atom " << iAt->id << " is no hybrid"
                        << std::endl;
            }
            
            if (iAt->chiralIdx==1)
            {
                std::cout << "Atom " << iAt->id << " may be in a positive chiral center "
                        << std::endl;
            }
            else if (iAt->chiralIdx==-1)
            {
                std::cout << "Atom " << iAt->id << " may be in a negative chiral center "
                        << std::endl;
            }
            else if (iAt->chiralIdx==2)
            {
                std::cout << "Atom " << iAt->id 
                        << " may be in a chiral center but the volume sign undefined "
                        << std::endl;
            }
            else if (iAt->chiralIdx==0)
            {
                std::cout << "Atom " << iAt->id 
                        << " is not a chiral center" << std::endl;
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
        setDefaultCoordGeos();
        
        ID metals[] = {"Li", "li", "Na", "na", "K",  "k",  "Rb", "rb", "Cs", "cs", "Fr", "fr",
                     "Be", "be", "Mg", "mg", "Ca", "ca", "Sr", "sr", "Ba", "ba", "Ra", "ra",
                     "Sc", "sc", "Y",  "y",
                     "B", "b", "Si", "si", "Ge", "ge", "As", "as", "Sb", "sb", "Te", "te", "Po", "po",
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
        
        MetalTable.assign(metals, metals+123);
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
                std::cout <<iA->id << " is a metal atom." << std::endl;
                
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
        
    }
    
    void DictCifFile::addMissHydroAtoms()
    {
        // check valence and covalent bonds for each atom, add missing atom if 
        // necessary
        PeriodicTable   aPTab; 
        
        std::vector<ID> allMetals;
        initMetalTab(allMetals);
        
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
        /*
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
        */
        // remember the residue name is ignored so we have
        // (int)tF.size() == (int)hasProps["angle"].size()+1
        //if ((int)tF.size() == (int)hasProps["angle"].size()+1 
        //        && (int)tF.size() >2 && tF[0].find("#") ==std::string::npos)
        //{
            // itsCurAngle = new AngleDict();
            // itsCurAngle->seriNum = itsCurAngleSeriNum;
            // itsCurAngleSeriNum++;
            
            
            /*
             * No need for that part, could rewrite them
            if (hasProps["angle"].find("resName") != hasProps["angle"].end())
            {
                itsCurAngle->resName=tF[hasProps["angle"]["resName"]];
            }
            
            
            if (hasProps["angle"].find("atom_id_1") != hasProps["angle"].end())
            {
                cleanSymbol(tF[hasProps["angle"]["atom_id_1"]], "\"");
                itsCurAngle->atoms.push_back(TrimSpaces(tF[hasProps["angle"]["atom_id_1"]]));
            }
            
            if (hasProps["angle"].find("atom_id_2") != hasProps["angle"].end())
            {
                cleanSymbol(tF[hasProps["angle"]["atom_id_2"]], "\"");
                itsCurAngle->atoms.push_back(TrimSpaces(tF[hasProps["angle"]["atom_id_2"]]));
            }
            
            if (hasProps["angle"].find("atom_id_3") != hasProps["angle"].end())
            {
                cleanSymbol(tF[hasProps["angle"]["atom_id_3"]], "\"");
                itsCurAngle->atoms.push_back(TrimSpaces(tF[hasProps["angle"]["atom_id_3"]]));
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
            */
            //allAngles.push_back(*itsCurAngle);
            //delete itsCurAngle;
            //itsCurAngle = NULL;
        //}   
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
            REAL iniValue;
            
            int n = 0;
            
            for (int i=0; i <(int)tV1.size(); i++)
            {
                if(tS1 !=-1 && tS2 !=-1)
                {
                    if (tF=="even")
                    {
                        iniValue= 60.0 +  n*perValue;
                    }
                    else if (tF=="old")
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
                                              RingDict & tR)
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
        
        std::string clibMonDir(std::getenv("CLIBD_MON"));
        std::string metDefCoordGeoFileName = clibMonDir + "/allMetalDefCoordGeos.table";
        //std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
        //std::string metDefCoordGeoFileName = clibMonDir + "/lib/allMetalDefCoordGeos.table";
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
    
    extern void outMMCif(FileName tFName, 
                         ID tMonoRootName,
                         std::vector<LIBMOL::AtomDict>& tAtoms,
                         std::vector<int>    & tHydroAtoms,
                         std::vector<LIBMOL::BondDict>& tBonds, 
                         std::vector<LIBMOL::AngleDict>& tAngs, 
                         std::vector<LIBMOL::TorsionDict>& tTorsions, 
                         std::map<LIBMOL::ID, std::vector<LIBMOL::RingDict> > & tRings, 
                         std::vector<LIBMOL::PlaneDict>& tPlas, 
                         std::vector<LIBMOL::ChiralDict>& tChs)
    {
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
                    << "_chem_comp.number_atoms_nh" << std::endl
                    << "_chem_comp.desc_level" << std::endl
                    << tMonoRootName <<"\t"<< tMonoRootName << "\t" << "'.\t\t'\t"
                    << "non-polymer\t" << (int)tAtoms.size() << "\t" 
                    << (int)tAtoms.size()-(int)tHydroAtoms.size() << "\t."
                    << std::endl;
            
            outRestrF <<"# ------------------------------------------------------" << std::endl
                      <<"# ------------------------------------------------------" << std::endl
                      <<"#" << std::endl
                      <<"# --- DESCRIPTION OF MONOMERS ---" << std::endl
                      <<"#" << std::endl
                      <<"data_comp_" << tMonoRootName << std::endl
                      <<"#" << std::endl; 
        
            if (tAtoms.size() >0)
            {
                // atom info section           
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
                    outRestrF << tMonoRootName  
                              << std::setw(6) << iA->id 
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
            
            if (tBonds.size() >0)
            {
                // Bond sections 
                outRestrF << "loop_" << std::endl
                          << "_chem_comp_bond.comp_id" << std::endl
                          << "_chem_comp_bond.atom_id_1" << std::endl
                          << "_chem_comp_bond.atom_id_2" << std::endl
                          << "_chem_comp_bond.type" << std::endl
                          << "_chem_comp_bond.value_dist"<< std::endl
                          << "_chem_comp_bond.value_dist_esd" << std::endl;
                          // << "_chem_comp_bond.exact_cod_dist" << std::endl;
                
                for (std::vector<BondDict>::iterator iB=tBonds.begin();
                          iB !=tBonds.end(); iB++)
                {
                    outRestrF << tMonoRootName  
                              << std::setw(6)  << iB->atoms[0]  
                              << std::setw(6)  << iB->atoms[1]  
                              << std::setw(10) << iB->order 
                              << std::setw(10) << std::setprecision(3)
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
                    
                    if (tAtoms[iA->atoms[0]].isMetal)
                    {
                        for (std::vector<REAL>::iterator iCA=iA->codAngleValues.begin();
                                iCA !=iA->codAngleValues.end(); iCA++)
                        {
                            outRestrF << tMonoRootName << std::setw(6)
                                      << tAtoms[iA->atoms[1]].id << std::setw(6)
                                      << tAtoms[iA->atoms[0]].id << std::setw(6)
                                      << tAtoms[iA->atoms[2]].id ;
                            outRestrF << std::setw(6) << std::setprecision(3) <<  *iCA << "    "
                                      << std::setw(6) << std::setprecision(2) << iA->sigValue << std::endl;
                        }
                    }
                    else
                    {
                        outRestrF << tMonoRootName 
                                  << std::setw(6) << tAtoms[iA->atoms[1]].id
                                  << std::setw(6) << tAtoms[iA->atoms[0]].id 
                                  << std::setw(8) << tAtoms[iA->atoms[2]].id;
                        outRestrF << std::setw(10) << std::setprecision(3) <<  iA->value
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
                
                    outRestrF << tMonoRootName 
                              << std::setw(16) << iT->id
                              << std::setw(8)  << tAtoms[iT->atoms[0]].id 
                              << std::setw(8)  << tAtoms[iT->atoms[1]].id 
                              << std::setw(8)  << tAtoms[iT->atoms[2]].id 
                              << std::setw(8)  << tAtoms[iT->atoms[3]].id 
                              << std::setw(12) << std::setprecision(3) << iT->value  
                              << std::setw(8)  << "10.00" 
                              << std::setw(6)  << iT->period << std::endl;
                    idxTor++;        
                }
                
            }
            
            //  For chiral centers
            bool l_ch = false;
            if ((int)tChs.size() !=0)
            {
                l_ch = true;
            }
            else
            {
                for (int i_ch =0; i_ch < (int)tAtoms.size(); i_ch++)
                {
                    if (tAtoms[i_ch].chiralIdx ==1)
                    {
                        l_ch = true;
                        break;
                    }
                }
            }
            
            if (l_ch)
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
                    outRestrF << tMonoRootName << "    " 
                              << iCh->id  << "    ";
                    for (std::vector<int>::iterator iAt=iCh->atoms.begin();
                            iAt != iCh->atoms.end(); iAt++)
                    {
                        outRestrF << tAtoms[*iAt].id << "    ";
                    }
                    outRestrF << iCh->sign << std::endl;
                }
                // New chiral that are not in the input list 
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
                                  
                                outRestrF << tMonoRootName 
                                          << std::setw(10) << idxStr 
                                          << std::setw(6)  << iA->id;
                       
                                outRestrF << std::setw(6) << chirAtms[0] << "    ";
                                outRestrF << std::setw(6) << chirAtms[1] << "    ";
                                outRestrF << std::setw(6) << chirAtms[2] << "    ";
                                outRestrF << std::setw(6) << "BOTH" << std::endl;
                                idxC++;
                            }
                        }
                    }
                }
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
                            outRestrF << tMonoRootName 
                                      << std::setw(10) << idxPStr
                                      << std::setw(6)  << iAt->first 
                                      << std::setw(8)  << "0.020" << std::endl;
                        }
                        idxP++;
                    }
                }
            }
            outRestrF.close();
        }
    }
    
    
}
