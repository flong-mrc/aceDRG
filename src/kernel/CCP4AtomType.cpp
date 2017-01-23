/* 
 * File:   CCP4AtomType.cpp
 * Author: flong
 *
 * Created on September 6, 2012, 11:32 AM
 * 
 */

#include "CCP4AtomType.h"

namespace LIBMOL
{
    /* CPP4 atom energy types are built after atomic connections, whether they are in rings, 
     * and whether they are metal elements, have been determined.
     */
    
    CCP4AtomType::CCP4AtomType()
    {
    }
    
    CCP4AtomType::CCP4AtomType(const CCP4AtomType & tC)
    {     
    }
    
    CCP4AtomType::CCP4AtomType(const std::vector<AtomDict>& tAllAtoms, 
                               const std::map<ID, std::vector<RingDict> > & tAllRings)
    {
        for (std::vector<AtomDict>::const_iterator iA=tAllAtoms.begin();
                iA != tAllAtoms.end(); iA++)
        {
            allAtoms.push_back(*iA);
        }
        
        for (std::map<ID, std::vector<RingDict> >::const_iterator iR=tAllRings.begin();
                iR !=tAllRings.end(); iR++)
        {
            for (std::vector<RingDict>::const_iterator iRi=iR->second.begin();
                    iRi !=iR->second.end(); iRi++)
            {
                allRings[iR->first].push_back(*iRi);
            }
        }
    }
    
    CCP4AtomType::CCP4AtomType(const std::vector<AtomDict>& tAllAtoms, 
                               const std::vector<RingDict>& tAllRingsV)
    {
        for (std::vector<AtomDict>::const_iterator iA=tAllAtoms.begin();
                iA != tAllAtoms.end(); iA++)
        {
            allAtoms.push_back(*iA);
        }
        
        
        for (std::vector<RingDict>::const_iterator iR=tAllRingsV.begin();
                 iR !=tAllRingsV.end(); iR++)
        {
            allRingsV.push_back(*iR);
        }
        
    }
    
    CCP4AtomType::~CCP4AtomType()
    {
    }
    
    void CCP4AtomType::setAllAtomsCCP4Type()
    {
        PeriodicTable aPeriTab;
        
        // Three rounds (1) basic, (2) ring and chiral modify, 
        // (3) neighbor atom combination. 
        
        
        // CCP4 types for H atoms are not setup in the first round
        // Other atom will be considered with ring and chiral information
        
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA !=allAtoms.end(); iA++)
        {
            // Give the default ccp4 types to be element symbols
            if (iA->ccp4Type==NullString)
            {
                iA->ccp4Type=iA->chemType;
            }
            
            // set atom ccp4 types to be individual one
            if (aPeriTab.elements[iA->chemType]["matType"] !=1)
            {
                setOneAtomCCP4Type(aPeriTab, *iA);
            }
        }
        
        
        
        // Combining all neighboring atoms information
        // But do not do anything on metal atoms at present.
        /*
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA !=allAtoms.end(); iA++)
        {
            if (aPeriTab.elements[iA->chemType]["matType"] ==2 
                    || aPeriTab.elements[iA->chemType]["matType"]==8)
            {
                setOneAtomCCP4Type(aPeriTab, *iA);
            }
        }
         */
        // Now Hydrogen atom types
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA !=allAtoms.end(); iA++)
        {
            if (aPeriTab.elements[iA->chemType]["matType"] ==1) 
            {
                setOneAtomCCP4Type(aPeriTab, *iA);
            }
        }
    }
    
    void CCP4AtomType::setOneAtomCCP4Type(PeriodicTable & tP, 
                                          AtomDict& tAtom)
    {
        int nType = tP.elements[tAtom.chemType]["matType"];
        
        switch(nType)
        {
            case 1:
                setHydroAtomCCP4Type(tAtom);
                break;
            case 2:
                setOrgAtomCCP4Type(tAtom);
                break;
            case 3:
                SetAlkaliMetalsAtomCCP4Type(tAtom);
                break;
            case 4:
                SetAlkalineEarthMetalsAtomCCP4Type(tAtom);
                break;
            case 5:
                SetTransitionMetalsAtomCCP4Type(tAtom);
                break;
            case 6:
                SetOtherMetalAtomCCP4Type(tAtom);
                break;
            case 7:
                SetSemimetallicsAtomCCP4Type(tAtom);
                break;
            case 8:
                SetHalogensAtomCCP4Type(tAtom);
                break;
            case 9:
                SetRareEarthAtomCCP4Type(tAtom);
                break;
            case 10:
                SetInertGasesAtomCCP4Type(tAtom);
                break;
            default:
                std::cout << "Error: Which material type is this element,  "
                        << "Atom ID: " << tAtom.id 
                        << " its material type (should be smaller than 10) : "
                        << nType << std::endl;
                exit(1);        
        }
        
        StrUpper(tAtom.ccp4Type);
    }
    
    void CCP4AtomType::setHydroAtomCCP4Type(AtomDict& tAtom)
    {
        // Hydrogen atom types are decided in the last step.
        // Only covalence bonds to H are considered for CCP4 type 
        //in this stage. Hydrogen bonds will be added in sometime later
        
        tAtom.ccp4Type = "H";
        
        if ((int)tAtom.connAtoms.size() == 1)   
        {
            int i = tAtom.connAtoms[0];
            //std::cout << "NB atom id " << allAtoms[i].id << std::endl;
            //std::cout << "NB atom element type " << allAtoms[i].chemType << std::endl;
            if (allAtoms[i].chemType.compare("S")==0)
            {
                tAtom.ccp4Type="HSH1";
            }
        }
        /*
        if ((int)tAtom.connAtoms.size() == 1 && (int)tAtom.ccp4Type.size()==0)   
        {
            int i = tAtom.connAtoms[0];
            
            if (allAtoms[i].chemType.compare("C")==0)
            {
                ID cNum=allAtoms[i].id.substr(1); // e.g. 11 for C11, 25 for C25 
                if (allAtoms[i].connHAtoms.size()==1)
                {
                    tAtom.ccp4Type = tAtom.chemType +  cNum; 
                }
                else
                {
                    // C allows max 6 NB Hs
                    std::vector<ID> letts;
                    letts.push_back("A"); 
                    letts.push_back("B");
                    letts.push_back("C");
                    letts.push_back("D");
                    letts.push_back("E");
                    letts.push_back("F");
                    
                    for (int nNB=0; nNB <(int)allAtoms[i].connHAtoms.size();
                            nNB++)
                    {
                        int j = allAtoms[i].connHAtoms[nNB];
                        allAtoms[j].ccp4Type = "H" + cNum + letts[nNB];
                    }
                }
            }
            tAtom.ccp4Type = tAtom.chemType + allAtoms[i].ccp4Type; 
        }
        else
        {
            int i = tAtom.connAtoms[0];
            if (allAtoms[i].connHAtoms.size()==1)
            {
                tAtom.ccp4Type = tAtom.chemType + allAtoms[i].id; 
            }
            else
            {
                for (int nNB=0; nNB <(int)allAtoms[i].connHAtoms.size();
                            nNB++)
                {
                    int j = allAtoms[i].connHAtoms[nNB];
                    allAtoms[j].ccp4Type = "H" + allAtoms[j].id + IntToStr(nNB);
                }
            }
        }
        */
        
        
    }
    
    void CCP4AtomType::setOrgAtomCCP4Type(AtomDict& tAtom)
    {
        // Residues of charged group
        std::map<ID, REAL> chargedGrp;
        chargedGrp["ARG"] = 0.0;   // value of charge
        
        int R5 = 0;
        int R6 = 0;
            
        for (std::map<ID, int>::iterator iR=tAtom.ringRep.begin();
                    iR!=tAtom.ringRep.end(); iR++)
        {
            if (iR->second ==5)
            {
                R5+=1;
            }
            if(iR->second==6)
            {
                R6+=1;
            }
        }
        /*
        std::cout << "atom " << tAtom.id << std::endl
                  << "it is ring 5 " << R5 << std::endl
                  << " it is ring 6 " << R6 << std::endl
                  << "its bond-idx " << tAtom.bondingIdx << std::endl
                  << "its connect to " << (int)tAtom.connHAtoms.size() 
                  << " H atoms " << std::endl;
        */
        // C atom
        if (tAtom.chemType.compare("C")==0)
        {
            //std::cout << "atom id " << tAtom.id << " C atom " << std::endl;
            //std::cout << " H atom connection " << tAtom.connHAtoms.size() << std::endl;
            if (tAtom.bondingIdx==2)
            {
                if (R5 && R6)
                {
                    tAtom.ccp4Type = "CR56";
                }
                else if (R5==2)
                {
                    tAtom.ccp4Type ="CR55";
                }
                else if (R6==2)
                {
                    tAtom.ccp4Type = "CR66";
                }
                else if(R5==1)
                {
                    if((int)tAtom.connHAtoms.size()==1)
                    {
                        tAtom.ccp4Type ="CR15";
                    }
                    else if((int)tAtom.connHAtoms.size()==0)
                    {
                        tAtom.ccp4Type ="CR5";
                    }
                }
                else if(R6==1)
                {
                    if((int)tAtom.connHAtoms.size()==1)
                    {
                        tAtom.ccp4Type ="CR16";
                    }
                    else if((int)tAtom.connHAtoms.size()==0)
                    {
                        tAtom.ccp4Type ="CR6";
                    }
                }
                else  // No Rings at all
                {
                    if ((int)tAtom.connHAtoms.size() == 1)
                    {
                        tAtom.ccp4Type = "C1";
                    }
                    else if ((int)tAtom.connHAtoms.size() == 2)
                    {
                        tAtom.ccp4Type = "C2";
                    }
                    else if ((int)tAtom.connHAtoms.size() == 0)
                    {
                        tAtom.ccp4Type = "C";
                    }
                }
            }
            else if (tAtom.bondingIdx==3)
            {
                // SP3 bonding
                if ((int)tAtom.connHAtoms.size() ==0)
                {
                    tAtom.ccp4Type = "CT";
                }
                else if ((int)tAtom.connHAtoms.size() == 1)
                {
                    tAtom.ccp4Type = "CH1";
                }
                else if ((int)tAtom.connHAtoms.size() == 2)
                {
                    tAtom.ccp4Type = "CH2";
                }
                else if ((int)tAtom.connHAtoms.size() == 3)
                {
                    tAtom.ccp4Type = "CH3";
                } 
            }
            else if (tAtom.bondingIdx==1)
            {
                tAtom.ccp4Type = "CSP";
            }
        }
        else if (tAtom.chemType.compare("N")==0)
        {
            // PI bonding
            
            // SP2 bonding
            if (tAtom.bondingIdx==2)
            {
                std::map<ID, REAL>::iterator iFC;
                iFC = chargedGrp.find(tAtom.resName);
                
                if (tAtom.connAtoms.size() == 3 && R5)
                {
                    tAtom.ccp4Type = "NR5";
                }
                else if (tAtom.connAtoms.size() == 3 && R6)
                {
                    tAtom.ccp4Type = "NR6";
                }
                else if (tAtom.connHAtoms.size() == 1)
                {
                    if (R5)
                    {
                        tAtom.ccp4Type = "NR15";
                    }
                    else if (R6)
                    {
                        tAtom.ccp4Type = "NR16";
                    }
                    else if (iFC !=chargedGrp.end())
                    {
                        tAtom.ccp4Type ="NC1";
                    }
                    else
                    {
                        tAtom.ccp4Type ="NH1";
                    }
                }
                else if (tAtom.connHAtoms.size() == 2)
                {
                    if (iFC !=chargedGrp.end())
                    {
                        tAtom.ccp4Type ="NC2";
                    }
                    else
                    {
                        tAtom.ccp4Type ="NH2";
                    }
                }
                else if (tAtom.connHAtoms.size() ==0)
                {
                    if (tAtom.connAtoms.size() == 2 && R6)
                    {
                        tAtom.ccp4Type = "NRD6";
                    }
                    else if (tAtom.connAtoms.size() == 2 && R5)
                    {
                        tAtom.ccp4Type = "NRD5";
                    }
                    else
                    {
                        tAtom.ccp4Type = "N";
                    }
                }
                
            }
            else if (tAtom.bondingIdx ==3)  // SP3
            {
                if (tAtom.connHAtoms.size()==1)
                {
                    tAtom.ccp4Type="NT1";
                }
                else if(tAtom.connHAtoms.size()==2)
                {
                    tAtom.ccp4Type = "NT2";
                }
                else if(tAtom.connHAtoms.size()==3)
                {
                    tAtom.ccp4Type = "NT3";
                }
                else
                {
                    tAtom.ccp4Type = "NT";
                }
            }
            else if (tAtom.bondingIdx==1)
            {
                tAtom.ccp4Type = "NS";
            }
            
        }  
        else if (tAtom.chemType.compare("P")==0)
        {
            if((int)tAtom.connAtoms.size() ==4)
            {
                tAtom.ccp4Type="P";
            }
            else
            {
                tAtom.ccp4Type ="P1";
            }
        }
        else if (tAtom.chemType.compare("O")==0)
        {
            bool lP=false, lS=false, lB=false;
            for (std::vector<int>::iterator iNB=tAtom.connAtoms.begin();
                             iNB!=tAtom.connAtoms.end(); iNB++)
            {
                if (allAtoms[*iNB].chemType.compare("P")==0)
                {
                    lP=true;
                }
                if (allAtoms[*iNB].chemType.compare("S")==0)
                {
                    lS=true;
                }
                if (allAtoms[*iNB].chemType.compare("B")==0)
                {
                    lB=true;
                }
            }
            
            // SP2 bonding
            if (tAtom.bondingIdx==2)
            {
                // SP2 bonding
                //std::cout << "O atom " << " charge " << tAtom.parCharge << std::endl
                //          << "connections " << tAtom.connAtoms.size() << std::endl;
                if (tAtom.parCharge)
                { 
                    if(lP)
                    {
                        tAtom.ccp4Type = "OP";
                    }
                    else if(lS)
                    {
                        tAtom.ccp4Type = "OS";
                    }
                    else if(lB)
                    {
                        tAtom.ccp4Type = "OB";
                    }
                    else
                    {
                        tAtom.ccp4Type = "OC";
                    }
                }
                else if (tAtom.connAtoms.size() ==2)
                {
                    if ((int)tAtom.connHAtoms.size()==1)
                    {
                        // O in alcohol groups
                        tAtom.ccp4Type = "OH1";
                    }
                    else if ((int)tAtom.connHAtoms.size()==2)
                    {
                        // O in for example water molecules 
                        tAtom.ccp4Type = "OH2";
                    }
                    else
                    {
                        tAtom.ccp4Type = "O2";
                    }
                }  
                else
                {
                    if (tAtom.formalCharge)
                    {
                        if(lP)
                        {
                            tAtom.ccp4Type = "OP";
                        }
                        else if(lS)
                        {
                            tAtom.ccp4Type = "OS";
                        }
                        else if(lB)
                        {
                            tAtom.ccp4Type = "OB";
                        }
                        else
                        {
                            tAtom.ccp4Type = "OC";
                        }
                    }
                    else
                    {
                        tAtom.ccp4Type = "O";
                    }
                }
            }
            else if (tAtom.bondingIdx==3)
            {
                bool lC=false;
                for (std::vector<int>::iterator iNB=tAtom.connAtoms.begin();
                                iNB!=tAtom.connAtoms.end(); iNB++)
                {
                    if (allAtoms[*iNB].chemType.compare("C")==0)
                    {
                        lC=true;
                    }
                }
                if (lC && (int)tAtom.connHAtoms.size()==1 
                       && (int)tAtom.connAtoms.size()==2)
                {
                    // O in alcohol groups
                    tAtom.ccp4Type = "OH1";
                }
                else if (tAtom.connHAtoms.size()==2)
                {
                    // oxygen of water
                    tAtom.ccp4Type = "OH2";
                }
                // oxygen of water in MO6 ?
                else if (tAtom.connAtoms.size() ==2)
                {
                    if(tAtom.parCharge)
                    {
                        tAtom.ccp4Type = "OC2";
                    }
                    else if (tAtom.connHAtoms.size()==1)
                    {
                        // O in alcohol groups
                        tAtom.ccp4Type = "OH1";
                    }
                    else
                    {
                        tAtom.ccp4Type = "O2";
                    }
                }  
            }
            else if ((int)tAtom.connAtoms.size() ==1)
            {
                /*
                bool lP=false, lS=false, lB=false;
                for (std::vector<int>::iterator iNB=tAtom.connAtoms.begin();
                             iNB!=tAtom.connAtoms.end(); iNB++)
                {
                    if (allAtoms[*iNB].chemType.compare("P")==0)
                    {
                        lP=true;
                    }
                    if (allAtoms[*iNB].chemType.compare("S")==0)
                    {
                            lS=true;
                    }
                    if (allAtoms[*iNB].chemType.compare("B")==0)
                    {
                        lB=true;
                    }
                    
                }
                        
                if(lP)
                {
                    tAtom.ccp4Type = "OP";
                }
                else if(lS)
                {
                    tAtom.ccp4Type = "OS";
                }
                else if(lB)
                {
                    tAtom.ccp4Type = "OB";
                }
                else
                {
                    tAtom.ccp4Type = "OC";
                }
                 */
                if (tAtom.formalCharge)
                    {
                        if(lP)
                        {
                            tAtom.ccp4Type = "OP";
                        }
                        else if(lS)
                        {
                            tAtom.ccp4Type = "OS";
                        }
                        else if(lB)
                        {
                            tAtom.ccp4Type = "OB";
                        }
                        else
                        {
                            tAtom.ccp4Type = "OC";
                        }
                    }
                    else
                    {
                        tAtom.ccp4Type = "O";
                    }
            }
            
            else
            {
                tAtom.ccp4Type = "O";
            }
        }
        else if (tAtom.chemType.compare("S")==0)
        {
            if (tAtom.connAtoms.size() ==3 || tAtom.connAtoms.size() ==4)
            {
                if(tAtom.connHAtoms.size()==0)
                {
                    tAtom.ccp4Type = "S3";
                }
                else if (tAtom.connHAtoms.size()==1)
                {
                    tAtom.ccp4Type = "SH1";
                }
            }
            else if (tAtom.connAtoms.size() ==2)
            {
                if(tAtom.connHAtoms.size()==0)
                {
                    tAtom.ccp4Type = "S2";
                }
                else
                {
                    tAtom.ccp4Type = "SH1";
                }
                
                // should not have the case that 2 connections, both of them H atoms
            }
            else if ((int)tAtom.connAtoms.size() ==1)
            {
                // One double bond
                tAtom.ccp4Type = "S1";
            }
            else if(tAtom.connHAtoms.size()==0)
            {
                tAtom.ccp4Type = "S";
            }
            else if(tAtom.connHAtoms.size()==1)
            {
                tAtom.ccp4Type = "SH1";
            }
            else
            {
                tAtom.ccp4Type = "S";
            }
            // 2 H connections or above ?
        }
        else
        {
            // elements have no special treatment at the moment 
            tAtom.ccp4Type = tAtom.chemType;   
            StrUpper(tAtom.ccp4Type);
        }
        
        
        //std::cout << "Atom Name " << tAtom.id << " bonding index " 
        //              << tAtom.bondingIdx << std::endl;
        //std::cout << "Its element type " << tAtom.chemType << std::endl;
        
        //if (!tAtom.codClass.empty())
        //{
            //std::cout << " Cod type" << tAtom.codClass << std::endl; 
        //}
        //else
        //{
            //std::cout << std::endl;
        //}
        //std::cout << " Its ccp4 type " << tAtom.ccp4Type << std::endl;
    }
    
    void CCP4AtomType::SetAlkaliMetalsAtomCCP4Type(AtomDict& tAtom)
    {
        // Currently atoms of metal elements just take their chemType as CCP4 type
        // Should do more later on
        tAtom.ccp4Type = tAtom.chemType;
    }
    
    void CCP4AtomType::SetAlkalineEarthMetalsAtomCCP4Type(AtomDict& tAtom)
    {
        // Currently atoms of metal elements just take their chemType as CCP4 type
        // Should do more later on
        tAtom.ccp4Type = tAtom.chemType;
    }
    
    void CCP4AtomType::SetTransitionMetalsAtomCCP4Type(AtomDict& tAtom)
    {
        // Currently atoms of metal elements just take their chemType as CCP4 type
        // Should do more later on
        tAtom.ccp4Type = tAtom.chemType;
    }
    
    void CCP4AtomType::SetOtherMetalAtomCCP4Type(AtomDict& tAtom)
    {
        // Currently atoms of metal elements just take their chemType as CCP4 type
        // Should do more later on
        tAtom.ccp4Type = tAtom.chemType;
    }
    
    void CCP4AtomType::SetSemimetallicsAtomCCP4Type(AtomDict& tAtom)
    {
        // Currently atoms of metal elements just take their chemType as CCP4 type
        // Should do more later on
        tAtom.ccp4Type = tAtom.chemType;
    }
    
    void CCP4AtomType::SetHalogensAtomCCP4Type(AtomDict& tAtom)
    {
        // Currently atoms of metal elements just take their chemType as CCP4 type
        // Should do more later on
        tAtom.ccp4Type = tAtom.chemType;
    }
    
    void CCP4AtomType::SetRareEarthAtomCCP4Type(AtomDict& tAtom)
    {
        // Currently atoms of Rare-Earth elements just take their chemType as CCP4 type
        // Should do more later on
        tAtom.ccp4Type = tAtom.chemType;
    }
    
    void CCP4AtomType::SetInertGasesAtomCCP4Type(AtomDict& tAtom)
    {
        // Currently atoms of Inert-Gases elements just take their chemType as CCP4 type
        // Should do more later on
        tAtom.ccp4Type = tAtom.chemType;
    }
    
    CCP4DictParas::CCP4DictParas()
    {
        std::string clibMonDir(std::getenv("CLIBD_MON"));
        std::string fName(clibMonDir);
        fName.append("ener_lib.cif");
        
        std::ifstream fParams(fName.c_str());
        if (fParams.is_open())
        {
            bool iStart =false;
            std::string tRecord="";
            std::vector<std::string>  tBLs;
            std::vector<std::vector<std::string> > tBs;
            while (!fParams.eof())
            {
                std::getline(fParams, tRecord);
                
                if (tRecord.find("loop_") !=std::string::npos)
                {
                    if (tBLs.size() !=0)
                    {
                        tBs.push_back(tBLs);
                    }
                    tBLs.clear();
                }
                else if (iStart && (int)tRecord.size() !=0 && tRecord[0] !='#')
                {
                    tBLs.push_back(tRecord);
                }
                else if (tRecord.find("data_energy") !=std::string::npos)
                {
                    iStart = true;
                }
            }
            if (tBLs.size() !=0)
            {
                tBs.push_back(tBLs);
            }
            
            fParams.close();
            
            if (tBs.size() !=0)
            {
                for (std::vector<std::vector<std::string> >::iterator iBLs=tBs.begin();
                        iBLs != tBs.end(); iBLs++)
                {
                    for (std::vector<std::string>::iterator iL=iBLs->begin();
                            iL!=iBLs->end(); iL++)
                    {
                        // std::cout << *iL << std::endl;
                   
                        if (iL->find("_lib_atom.") != std::string::npos)
                        {   
                            getAtomPropsTable(iBLs);
                            break;
                        }
                        else if (iL->find("_lib_bond.") != std::string::npos)
                        {
                            getBondPropsTable(iBLs);
                            break;
                        }
                        else if (iL->find("_lib_angle.") != std::string::npos)
                        {
                            
                            getAnglePropsTable(iBLs);
                            break;
                        }
                    }
                }
                
            }
        }
        else
        {
            std::cout << fName << " can not be open for reading " << std::endl;
            exit(1);
        }
    }
    
    CCP4DictParas::~CCP4DictParas()
    {
    }
    
    void CCP4DictParas::getAtomPropsTable(std::vector<std::vector<std::string> >::iterator tLines)
    {
        // value -1 means parameter value not available
        REAL a=-1;
        std::map<std::string, int> propIdxs;
        int iLine =0;
        for (std::vector<std::string>::iterator iLi=tLines->begin();
                iLi !=tLines->end(); iLi++)
        {
            std::vector<std::string>            tBuf;
            StrTokenize(TrimSpaces(*iLi), tBuf);
            // std::cout << *iLi << std::endl; 
            if (tBuf.size()==1 && tBuf[0] !="#")
            {
                std::vector<std::string>      tBuf2;
                StrTokenize(TrimSpaces(*iLi), tBuf2, '.');
                if (tBuf2.size()==2)
                {
                    propIdxs[tBuf2[1]]=iLine;
                    iLine++;
                }
            }
            else if (tBuf.size() >= 2 && tBuf.size() == propIdxs.size() 
                     && tBuf[0] != ".")
            {
                for (std::map<std::string, int>::iterator iP=propIdxs.begin();
                        iP != propIdxs.end(); iP++)
                {
                    if (iP->first !="type"&& iP->first !="element")
                    {
                        if (tBuf[iP->second] =="." || tBuf[iP->second] =="N")     
                        {
                            atomPropsTable[tBuf[0]][iP->first] = a;
                        }
                        else if (tBuf[iP->second]=="D")
                        {
                            atomPropsTable[tBuf[0]][iP->first] = 1.0;
                        }
                        else if (tBuf[iP->second]=="A")
                        {
                            atomPropsTable[tBuf[0]][iP->first] = 2.0;
                        }
                        else if (tBuf[iP->second]=="B")
                        {
                            atomPropsTable[tBuf[0]][iP->first] = 3.0;
                        }
                        else if (tBuf[iP->second]=="H")
                        {
                            atomPropsTable[tBuf[0]][iP->first] = 4.0;
                        }
                        else
                        {
                            atomPropsTable[tBuf[0]][iP->first] = StrToReal(tBuf[iP->second]);
                        }
                    }
                    if (iP->first=="element")
                    {
                        atomTypeElementTable[tBuf[0]] = tBuf[iP->second];
                    }
                }
            }  
        }
        /*
        std::cout << "The following properties are associated with CCP4 atom types " << std::endl;
        for (std::map<std::string, std::map<std::string, REAL> >::iterator iAts=atomPropsTable.begin();
                iAts != atomPropsTable.end(); iAts++)
        {
            std::cout << "For atom type " << iAts->first << " : " << std::endl;
            std::cout << "Element = " << atomTypeElementTable[iAts->first] << std::endl;
            for ( std::map<std::string, REAL>::iterator iProp=iAts->second.begin();
                    iProp != iAts->second.end(); iProp++)
            {
                std::cout << iProp->first << " = " << iProp->second << std::endl;
            }
            std::cout << std::endl;
        }
        */        
        
    }
    
    void CCP4DictParas::getBondPropsTable(std::vector<std::vector<std::string> >::iterator tLines)
    {
        // value -1 means parameter value not available
        REAL a =-1.0;
        
        std::map<std::string, int> propIdxs;
        int iLine =0;
        for (std::vector<std::string>::iterator iLi=tLines->begin();
                iLi !=tLines->end(); iLi++)
        {
            std::vector<std::string>            tBuf;
            StrTokenize(TrimSpaces(*iLi), tBuf);
            if (tBuf.size()==1 && tBuf[0] !="#")
            {
                std::vector<std::string>      tBuf2;
                StrTokenize(TrimSpaces(*iLi), tBuf2, '.');
                if (tBuf2.size()==2)
                {
                    propIdxs[tBuf2[1]]=iLine;
                    iLine++;
                }
            }
            else if (tBuf.size() >= 2 && tBuf.size() == propIdxs.size() 
                     && tBuf[0] != ".")
            {
                if (tBuf[propIdxs["type"]]=="single")
                {
                    a = 1.0;
                }
                else if (tBuf[propIdxs["type"]]=="double")
                {
                    a = 2.0;
                }
                else if (tBuf[propIdxs["type"]]=="triple")
                {
                    a = 3.0;
                }
                else if (tBuf[propIdxs["type"]]=="aromatic")
                {
                    a = 4.0;
                }
                else if (tBuf[propIdxs["type"]]=="aromat")
                {
                    a = 5.0;
                }
                else if (tBuf[propIdxs["type"]]=="deloc")
                {
                    a = 6.0;
                }
                else if (tBuf[propIdxs["type"]]=="metal")
                {
                    a = 7.0;
                }
                bondPropsTable[tBuf[propIdxs["atom_type_1"]]][tBuf[propIdxs["atom_type_2"]]]
                       ["type"] = a;
                
                
                a = -1.0;
                if (tBuf[propIdxs["const"]] != ".")
                {
                    bondPropsTable[tBuf[propIdxs["atom_type_1"]]][tBuf[propIdxs["atom_type_2"]]]
                                  ["const"] = StrToReal(tBuf[propIdxs["const"]]);
                }
                else
                {
                    bondPropsTable[tBuf[propIdxs["atom_type_1"]]][tBuf[propIdxs["atom_type_2"]]]
                                  ["const"] = a;
                }
                
                if (tBuf[propIdxs["length"]] != ".")
                {
                    bondPropsTable[tBuf[propIdxs["atom_type_1"]]][tBuf[propIdxs["atom_type_2"]]]
                                  ["length"] = StrToReal(tBuf[propIdxs["length"]]);
                }
                else
                {
                    bondPropsTable[tBuf[propIdxs["atom_type_1"]]][tBuf[propIdxs["atom_type_2"]]]
                                  ["length"] = a;
                }
                
                if (tBuf[propIdxs["value_esd"]] != ".")
                {
                    bondPropsTable[tBuf[propIdxs["atom_type_1"]]][tBuf[propIdxs["atom_type_2"]]]
                                  ["value_esd"] = StrToReal(tBuf[propIdxs["value_esd"]]);
                }
                else
                {
                    bondPropsTable[tBuf[propIdxs["atom_type_1"]]][tBuf[propIdxs["atom_type_2"]]]
                                  ["value_esd"] = a;
                }
                
            }
        }
        
        
        
    }
    
    void CCP4DictParas::getAnglePropsTable(std::vector<std::vector<std::string> >::iterator tLines)
    {
        // value -1 means parameter value not available
        REAL a =-1.0;
        
        std::map<std::string, int> propIdxs;
        int iLine =0;
        for (std::vector<std::string>::iterator iLi=tLines->begin();
                iLi !=tLines->end(); iLi++)
        {
            std::vector<std::string>            tBuf;
            StrTokenize(TrimSpaces(*iLi), tBuf);
            if (tBuf.size()==1 && tBuf[0] !="#")
            {
                std::vector<std::string>      tBuf2;
                StrTokenize(TrimSpaces(*iLi), tBuf2, '.');
                if (tBuf2.size()==2)
                {
                    propIdxs[tBuf2[1]]=iLine;
                    iLine++;
                }
            }
            else if (tBuf.size() >= 2 && tBuf.size() == propIdxs.size() 
                     && tBuf[0] != ".")
            {
                if (tBuf[propIdxs["const"]] != ".")
                {
                    anglePropsTable[tBuf[propIdxs["atom_type_1"]]][tBuf[propIdxs["atom_type_2"]]]
                     [tBuf[propIdxs["atom_type_3"]]]["const"] = StrToReal(tBuf[propIdxs["const"]]);
                }
                else
                {
                    anglePropsTable[tBuf[propIdxs["atom_type_1"]]][tBuf[propIdxs["atom_type_2"]]]
                     [tBuf[propIdxs["atom_type_3"]]]["const"] = a;
                }
                
                if (tBuf[propIdxs["const"]] != ".")
                {
                    anglePropsTable[tBuf[propIdxs["atom_type_1"]]][tBuf[propIdxs["atom_type_2"]]]
                     [tBuf[propIdxs["atom_type_3"]]]["value"] = StrToReal(tBuf[propIdxs["value"]]);
                }
                else
                {
                    anglePropsTable[tBuf[propIdxs["atom_type_1"]]][tBuf[propIdxs["atom_type_2"]]]
                     [tBuf[propIdxs["atom_type_3"]]]["value"] = a;
                }
            }
        }
    }
    
    
        
    AtomTypeTool::AtomTypeTool()
    {
    }
    
    AtomTypeTool::AtomTypeTool(FileName tFname, FileType tFType)
    {
        if (tFType==CIF)
        {
            
        }
        
    }
    
    AtomTypeTool::AtomTypeTool(std::vector<AtomDict>& tAtoms, 
                               std::vector<BondDict>& tBonds, 
                               std::map<ID,std::vector<RingDict> > & tRings)
    {
        for (std::vector<AtomDict>::const_iterator iAt=tAtoms.begin();
                iAt !=  tAtoms.end(); iAt++)
        {
            allAtoms.push_back(*iAt);
        }
        
        for (std::vector<BondDict>::const_iterator iBo=tBonds.begin();
                iBo !=tBonds.end(); iBo++)
        {
            allBonds.push_back(*iBo);
        }
        
        for (std::map<ID, std::vector<RingDict> > ::const_iterator iRS=tRings.begin();
                iRS !=tRings.end(); iRS++)
        {
            for (std::vector<RingDict>::const_iterator iR=iRS->second.begin();
                    iR !=iRS->second.end(); iR++)
            {
                allRings[iRS->first].push_back(*iR);
            }
        }
        
    }
    
    AtomTypeTool::~AtomTypeTool()
    {
    }
    
    
}
