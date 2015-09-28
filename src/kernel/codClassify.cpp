/* 
 * File:   codClassify.cpp
 * Author: flong
 *
 * Created on April 6, 2012, 5:48 PM
 */
#include "codClassify.h"
#include "PDBFile.h"
#include "TransCoord.h"

namespace LIBMOL
{
    
    CodClassify::CodClassify():wSize(1000),
                               libmolTabDir("")
    {    
        pPeriodictable = new PeriodicTable();
    }
    
    
    /*
    CodClassify::CodClassify(const CodClassify& tCodC)
    {
        for (std::vector<AtomDict>::const_iterator iA=tCodC.allAtoms.begin();
                iA != tCodC.allAtoms.end(); iA++)
        {
            allAtoms.push_back(*iA);
        }
        
        for (std::vector<BondDict>::const_iterator iB=tCodC.allBonds.begin();
                iB != tCodC.allBonds.end(); iB++)
        {
            allBonds.push_back(*iB);
        }
        
        for (std::vector<AngleDict>::const_iterator iAn = tCodC.allAngles.begin();
                iAn !=tCodC.allAngles.end(); iAn++)
        {
            allAngles.push_back(*iAn);
        }
        
        for (std::vector<TorsionDict>::const_iterator iTo = tCodC.allTorsions.begin();
                iTo != tCodC.allTorsions.end(); iTo++)
        {
            allTorsions.push_back(*iTo);
        }
        
        for (std::vector<PlaneDict>::const_iterator iPl=tCodC.allPlanes.begin();
                iPl != tCodC.allPlanes.end(); iPl++)
        {
            allPlanes.push_back(*iPl);
        }
        
        for (std::map<ID, std::vector<RingDict> >::const_iterator iRi=tCodC.allRings.begin();
                iRi != tCodC.allRings.end(); iRi++)
        {
            allRings[iRi->first].push_back(iRi->second[0]);
        }
    }
    */
    
    CodClassify::CodClassify(const std::vector<AtomDict>& tAtoms):wSize(1000), 
                                                                  libmolTabDir("")
                                                                  
    {
        
        pPeriodictable = new PeriodicTable();
        
        for (std::vector<AtomDict>::const_iterator iA=tAtoms.begin();
                iA != tAtoms.end(); iA++)
        {
            allAtoms.push_back(*iA);
        }
    }
 
    CodClassify::CodClassify(const DictCifFile & tCifObj,
                             std::string   tLibmolTabDir):wSize(1000)
                                                         
    {
        libmolTabDir = tLibmolTabDir;
        
        pPeriodictable = new PeriodicTable();
        
        for (std::vector<AtomDict>::const_iterator iA=tCifObj.allAtoms.begin();
                iA != tCifObj.allAtoms.end(); iA++)
        {
            allAtoms.push_back(*iA);
        }
        
        for (std::vector<int>::const_iterator iH = tCifObj.allHAtomIdx.begin();
                iH != tCifObj.allHAtomIdx.end(); iH++)
        {
            allHAtomIdx.push_back(*iH);
        }
        
        for (std::vector<BondDict>::const_iterator iB=tCifObj.allBonds.begin();
                iB != tCifObj.allBonds.end(); iB++)
        {
            allBonds.push_back(*iB);
        }
        
       
        /* Do not copy angles into the new object. will re-setup them when all
         * information on atoms and bonds has been generated 
         
        for (std::vector<AngleDict>::const_iterator iAng=tCifObj.allAngles.begin();
                iAng != tCifObj.allAngles.end(); iAng++)
        {
            allAngles.push_back(*iAng);
        }
        
        for (std::map<int, std::vector<std::vector<int> > >::const_iterator iAngIdx =
                tCifObj.allAnglesIdxs.begin(); 
                iAngIdx !=tCifObj.allAnglesIdxs.end(); iAngIdx++)
        {
            
        }
  
        
        for (std::vector<TorsionDict>::const_iterator iTor=tCifObj.allTorsions.begin();
                iTor != tCifObj.allTorsions.end(); iTor++)
        {
            allTorsions.push_back(*iTor);
        }
    
         */    
        
        for (std::vector<ChiralDict>::const_iterator iCh=tCifObj.allChirals.begin();
                iCh !=tCifObj.allChirals.end(); iCh++)
        {
            allChirals.push_back(*iCh);
        }
        
        for (std::vector<int>::const_iterator iH=tCifObj.allHydroAtoms.begin();
                iH !=tCifObj.allHydroAtoms.end(); iH++)
        {
            allHydroAtoms.push_back(*iH);
        }
        
        setupSystem();
        
                
    }
    
    CodClassify::CodClassify(const std::vector<Molecule>& tMols):wSize(1000),
                                                                 libmolTabDir("")
    {
        pPeriodictable = new PeriodicTable();
        
        for (std::vector<Molecule>::const_iterator iMol=tMols.begin();
                iMol !=tMols.end(); iMol++)
        {
            for (std::vector<AtomDict>::const_iterator iA=iMol->atoms.begin();
                    iA !=iMol->atoms.end(); iA++)
            {
                allAtoms.push_back(*iA);
            }
            for (std::vector<BondDict>::const_iterator iBo=iMol->bonds.begin();
                    iBo !=iMol->bonds.end(); iBo++)
            {
                allBonds.push_back(*iBo);
            }
        }
    }
    
    CodClassify::CodClassify(const AllSystem& tAllSys):wSize(1000),
                                                       libmolTabDir("")
    {
        pPeriodictable = new PeriodicTable();
        
        for (std::vector<AtomDict>::const_iterator iA=tAllSys.allAtoms.begin();
                iA != tAllSys.allAtoms.end(); iA++)
        {
            allAtoms.push_back(*iA);
        }
        
        for (std::vector<int>::const_iterator iH = tAllSys.allHAtomIdx.begin();
                iH != tAllSys.allHAtomIdx.end(); iH++)
        {
            allHAtomIdx.push_back(*iH);
        }
        
        for (std::vector<BondDict>::const_iterator iB=tAllSys.allBonds.begin();
                iB != tAllSys.allBonds.end(); iB++)
        {
            allBonds.push_back(*iB);
        }
        
       
        /* Do not copy angles into the new object. will re-setup them when all
         * information on atoms and bonds has been generated 
         
        for (std::vector<AngleDict>::const_iterator iAng=tCifObj.allAngles.begin();
                iAng != tCifObj.allAngles.end(); iAng++)
        {
            allAngles.push_back(*iAng);
        }
        
        for (std::map<int, std::vector<std::vector<int> > >::const_iterator iAngIdx =
                tCifObj.allAnglesIdxs.begin(); 
                iAngIdx !=tCifObj.allAnglesIdxs.end(); iAngIdx++)
        {
            
        }
        */
        
        for (std::vector<TorsionDict>::const_iterator iTor=tAllSys.allTorsions.begin();
                iTor != tAllSys.allTorsions.end(); iTor++)
        {
            allTorsions.push_back(*iTor);
        }
        
        
        setupSystem();
    }
    
    CodClassify::CodClassify(const std::vector<AtomDict>& tAtoms, 
                             const std::vector<int>& tHAtomIdx, 
                             const std::vector<BondDict>& tBonds, 
                             const std::vector<AngleDict>& tAngles, 
                             const std::vector<TorsionDict>& tTorsions, 
                             const std::vector<ChiralDict>& tChirals, 
                             const std::vector<PlaneDict>& tPlans, 
                             const std::map<ID,std::vector<RingDict> >& tRings,
                             std::string               tLibmolTabDir):wSize(1000) 
    {
        libmolTabDir = tLibmolTabDir; 
        pPeriodictable = new PeriodicTable();
        
        for (std::vector<AtomDict>::const_iterator iA=tAtoms.begin();
                iA != tAtoms.end(); iA++)
        {
            allAtoms.push_back(*iA);
        }
        
        for (std::vector<int>::const_iterator iH = tHAtomIdx.begin();
                iH != tHAtomIdx.end(); iH++)
        {
            allHAtomIdx.push_back(*iH);
        }
        
        for (std::vector<BondDict>::const_iterator iB=tBonds.begin();
                iB != tBonds.end(); iB++)
        {
            allBonds.push_back(*iB);
        }
        
        for (std::vector<AngleDict>::const_iterator iAng=tAngles.begin();
                iAng != tAngles.end(); iAng++)
        {
            allAngles.push_back(*iAng);
        }
        
        for (std::vector<TorsionDict>::const_iterator iTor=tTorsions.begin();
               iTor != tTorsions.end(); iTor++)
        {
            allTorsions.push_back(*iTor); 
        }
        
        for (std::vector<ChiralDict>::const_iterator iCh=tChirals.begin();
                iCh !=tChirals.end(); iCh++)
        {
            allChirals.push_back(*iCh);
        }
        
        for (std::vector<int>::const_iterator iH=tHAtomIdx.begin();
                iH !=tHAtomIdx.end(); iH++)
        {
            allHydroAtoms.push_back(*iH);
        }
        
        for (std::vector<PlaneDict>::const_iterator iP=tPlans.begin();
                iP!=tPlans.end(); iP++)
        {
            allPlanes.push_back(*iP);
        }
                
        for (std::map<ID, std::vector<RingDict> >::const_iterator iR=tRings.begin();
                iR !=tRings.end(); iR++)
        {
            for (std::vector<RingDict>::const_iterator iR1=iR->second.begin();
                    iR1 !=iR->second.end(); iR1++)
            {
                allRings[iR->first].push_back(*iR1);
            }
        }
        
        // setupSystem();
        setupSystem();
        
    }
        
    
    CodClassify::CodClassify(const std::vector<AtomDict>& tAtoms, 
                             const std::vector<int>& tHAtomIdx, 
                             const std::vector<BondDict>& tBonds, 
                             const std::vector<AngleDict>& tAngles, 
                             const std::vector<TorsionDict>& tTorsions, 
                             const std::vector<ChiralDict>& tChirals, 
                             const std::vector<PlaneDict>& tPlans, 
                             const std::map<ID,std::vector<RingDict> >& tRings,
                             std::string               tLibmolTabDir,
                             int                       nTM):wSize(1000) 
    {
        libmolTabDir = tLibmolTabDir; 
        pPeriodictable = new PeriodicTable();
        
        for (std::vector<AtomDict>::const_iterator iA=tAtoms.begin();
                iA != tAtoms.end(); iA++)
        {
            allAtoms.push_back(*iA);
        }
        
        for (std::vector<int>::const_iterator iH = tHAtomIdx.begin();
                iH != tHAtomIdx.end(); iH++)
        {
            allHAtomIdx.push_back(*iH);
        }
        
        for (std::vector<BondDict>::const_iterator iB=tBonds.begin();
                iB != tBonds.end(); iB++)
        {
            allBonds.push_back(*iB);
        }
        
        for (std::vector<AngleDict>::const_iterator iAng=tAngles.begin();
                iAng != tAngles.end(); iAng++)
        {
            allAngles.push_back(*iAng);
        }
        
        for (std::vector<TorsionDict>::const_iterator iTor=tTorsions.begin();
               iTor != tTorsions.end(); iTor++)
        {
            allTorsions.push_back(*iTor); 
        }
        
        for (std::vector<ChiralDict>::const_iterator iCh=tChirals.begin();
                iCh !=tChirals.end(); iCh++)
        {
            allChirals.push_back(*iCh);
        }
        
        for (std::vector<int>::const_iterator iH=tHAtomIdx.begin();
                iH !=tHAtomIdx.end(); iH++)
        {
            allHydroAtoms.push_back(*iH);
        }
        
        for (std::vector<PlaneDict>::const_iterator iP=tPlans.begin();
                iP!=tPlans.end(); iP++)
        {
            allPlanes.push_back(*iP);
        }
                
        for (std::map<ID, std::vector<RingDict> >::const_iterator iR=tRings.begin();
                iR !=tRings.end(); iR++)
        {
            for (std::vector<RingDict>::const_iterator iR1=iR->second.begin();
                    iR1 !=iR->second.end(); iR1++)
            {
                allRings[iR->first].push_back(*iR1);
            }
        }
        
        // setupSystem();
        setupSystem2();
        
    }
   
    
    CodClassify::~CodClassify()
    {
        if(pPeriodictable)
        {
            delete pPeriodictable;
            pPeriodictable = NULL;
        }
    }

    int CodClassify::atomPosition(ID tID)
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
    
    /*
    void CodClassify::setupSystem()
    {
        
        int cLev = 2;      // should be a constant from an input file
        
       
        codAtomClassify(cLev);
        
        setAtomsBondingAndChiralCenter();
        
        getCCP4BondAndAngles();
        
        // std::vector<RingDict>          tmpRings;
        // std::vector<std::vector<int> > tmpAtoms;
        
        if (allRings.size() !=0)
        {
            for (std::map<ID, std::vector<RingDict> >::iterator iMr=allRings.begin();
                    iMr !=allRings.end(); iMr++)
            {
                for(std::vector<RingDict>::iterator iR=iMr->second.begin();
                        iR !=iMr->second.end(); iR++)
                {
                    // tmpRings.push_back(*iR);
                    
                    iR->setPlaneProp();
                    
                    
                    std::string PL("No");
                    if (iR->isPlanar)
                    {
                        PL = "Yes";
                    }
                    
                    std::cout << "Is ring, " << iMr->first << ", a planar ring? "
                              <<  PL << std::endl;
                   
                }
            }
        }
        
        
        mergePlaneRings(tmpRings, tmpAtoms);
        
        std::cout << "Number of merged atom rings " 
                  << tmpAtoms.size() << std::endl;
        
        for (std::vector<std::vector<int> >::iterator iR=tmpAtoms.begin();
                iR !=tmpAtoms.end(); iR++)
        {
            std::cout << "Atoms in a merged ring : " << std::endl;
            for (std::vector<int>::iterator iAt=iR->begin();
                    iAt !=iR->end(); iAt++)
            {
                std::cout << allAtoms[*iAt].id << std::endl;
            }
        }
         
        exit(0);
        
        
        // detectPlaneGroups();
        
        // initTargetAngles();
    }
    */
    
    void CodClassify::setupSystem()
    {
        
        int cLev = 2;      // should be a constant from an input file
        
       
        codAtomClassify(cLev);
        
        // codAtomClassifyNew2(cLev);
        
        setAtomsBondingAndChiralCenter();
        
        getCCP4BondAndAngles();
        
        // std::vector<RingDict>          tmpRings;
        // std::vector<std::vector<int> > tmpAtoms;
        
        
        
        if (allRings.size() !=0)
        {
            for (std::map<ID, std::vector<RingDict> >::iterator iMr=allRings.begin();
                    iMr !=allRings.end(); iMr++)
            {
                for(std::vector<RingDict>::iterator iR=iMr->second.begin();
                        iR !=iMr->second.end(); iR++)
                {
                    // tmpRings.push_back(*iR);
                    
                    iR->setPlaneProp();
                    
                    //setSugarRingInitComf(allAtoms, allTorsions, iR);
                    
                    checkOneSugarRing(allAtoms, iR);
                    if (iR->isSugar.compare("pyranose")==0)
                    {    
                        std::cout << "Find one pyranose ring " << std::endl;
                        // A pyranose ring, set torsions within the ring    
                        setPyranoseChairComf(allAtoms, iR, allTorsions);
                    }
                   
                    /*
                    std::string PL("No");
                    if (iR->isPlanar)
                    {
                        PL = "Yes";
                    }
                    
                    std::cout << "Is ring, " << iMr->first << ", a planar ring? "
                              <<  PL << std::endl;
                    */
                }
            }
        }
        
        // Ring aromaticity has been set, put that and other in atom properties
        for (std::vector<AtomDict>::iterator iAt=allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {
            iAt->setBaseRingProps();
        }
        
        /*
        std::cout << "Torsion angles are : " << std::endl;
        for (std::vector<TorsionDict>::iterator iTor=allTorsions.begin();
                iTor !=allTorsions.end(); iTor++)
        {
            std::cout << "Formed by " << allAtoms[iTor->atoms[0]].id 
                      << "_" << allAtoms[iTor->atoms[1]].id  << "_"
                      << allAtoms[iTor->atoms[2]].id << "_" 
                      << allAtoms[iTor->atoms[3]].id << ", value is: "
                      << iTor->value << std::endl;
        }
        */
        
        
        
    }
    
    void CodClassify::setupSystem2()
    {
        
        int cLev = 2;      // should be a constant from an input file
       
        // codAtomClassify(cLev);
        
        codAtomClassifyNew2(cLev);
        
        setAtomsBondingAndChiralCenter();
        
        getCCP4BondAndAngles();
        
        // std::vector<RingDict>          tmpRings;
        // std::vector<std::vector<int> > tmpAtoms;
        
        if (allRings.size() !=0)
        {
            for (std::map<ID, std::vector<RingDict> >::iterator iMr=allRings.begin();
                    iMr !=allRings.end(); iMr++)
            {
                for(std::vector<RingDict>::iterator iR=iMr->second.begin();
                        iR !=iMr->second.end(); iR++)
                {
                    // tmpRings.push_back(*iR);
                    
                    iR->setPlaneProp();
                    
                    //setSugarRingInitComf(allAtoms, allTorsions, iR);
                    
                    checkOneSugarRing(allAtoms, iR);
                    if (iR->isSugar.compare("pyranose")==0)
                    {    
                        std::cout << "Find one pyranose ring " << std::endl;
                        // A pyranose ring, set torsions within the ring    
                        setPyranoseChairComf(allAtoms, iR, allTorsions);
                    }
                   
                    /*
                    std::string PL("No");
                    if (iR->isPlanar)
                    {
                        PL = "Yes";
                    }
                    
                    std::cout << "Is ring, " << iMr->first << ", a planar ring? "
                              <<  PL << std::endl;
                    */
                }
            }
        }
        
        // Ring aromaticity has been set, put that and other in atom properties
        for (std::vector<AtomDict>::iterator iAt=allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {
            iAt->setBaseRingProps();
            
        }
        
        
        for (std::vector<AtomDict>::iterator iAt=allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {
            std::cout << "Atom " << iAt->seriNum << " : " << std::endl
                      << "Its ID " << iAt->id << std::endl
                      << "Its element type " << iAt->chemType
                      << "Its acedrg atom type " << iAt->codClass << std::endl
                      << "Its acedrg atom main type " << iAt->codAtmMain << std::endl
                      << "Its ccp4 atom type " << iAt->ccp4Type 
                      << std::endl;
           
            
        }
        
        
        
        /*
        std::cout << "Torsion angles are : " << std::endl;
        for (std::vector<TorsionDict>::iterator iTor=allTorsions.begin();
                iTor !=allTorsions.end(); iTor++)
        {
            std::cout << "Formed by " << allAtoms[iTor->atoms[0]].id 
                      << "_" << allAtoms[iTor->atoms[1]].id  << "_"
                      << allAtoms[iTor->atoms[2]].id << "_" 
                      << allAtoms[iTor->atoms[3]].id << ", value is: "
                      << iTor->value << std::endl;
        }
        */
        
    }
    
    
    void CodClassify::codAtomClassify(int dLev)
    {
        
        ringDetecting();
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
      
        
        // std::cout <<std::endl << "Output Atom COD classes now " << std::endl << std::endl;
        
        for (int i=0; i < (int)allAtoms.size(); i++)
        {
            setAtomCodClassName2(allAtoms[i], allAtoms[i], dLev);
            std::cout <<std::endl << "For atom " << allAtoms[i].id << std::endl 
                      << "class is " << allAtoms[i].codClass << std::endl;            
        }
        
        // set a hashing code and primeNB symbol to each atom.
        
        hashingAtoms();
        
        setAtomsNBSymb();
        
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
    
    void CodClassify::codAtomClassify2(int dLev)
    {
        
        ringDetecting2();
       
        /*
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA !=allAtoms.end(); iA++)
        {
           
            if ((int)iA->ringRep.size() >0)
            {
                std::cout << "From ring rep view " << std::endl;
                std::cout << "Atom " << iA->id << " is in  " 
                        << (int)iA->ringRep.size() << " ring(s)" << std::endl;
                for (std::map<std::string, int>::iterator iRR=iA->ringRep.begin();
                        iRR != iA->ringRep.end(); iRR++)
                {
                    std::cout << "Ring: " << iRR->first << ", size : "
                            << iRR->second << std::endl;
                }
            }
            
            if (iA->inRings.size() >0)
            {
                std::cout << "From ring number view " << std::endl;
                std::cout << "Atom " << iA->id << " is in  " 
                          << (int)iA->inRings.size() << " ring(s)" << std::endl;
                
                for (std::vector<int>::iterator iI=iA->inRings.begin();
                         iI !=iA->inRings.end(); iI++)
                {
                    std::cout << "Ring: " << *iI << std::endl;
                }
            }
           
            if (iA->ringRep.size() !=iA->inRings.size())
            {
               exit(1);
            }
        }
       
        */
        
        
        
        // std::cout <<std::endl << "Output Atom COD classes now " << std::endl << std::endl;
        
        for (int i=0; i < (int)allAtoms.size(); i++)
        {
            setAtomCodClassName2(allAtoms[i], allAtoms[i], dLev);
            //std::cout <<std::endl << "For atom " << allAtoms[i].id << std::endl 
            //          << "class is " << allAtoms[i].codClass << std::endl;            
        }
        
        
        // set a hashing code and primeNB symbol to each atom.
        
        //hashingAtoms();
        
        //setAtomsNBSymb();
        
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
    
    void CodClassify::codAtomClassifyNew(int dLev)
    {
        ringTools aRingTool;
        allRings.clear();
        allPlanes.clear();
        int nMaxRing = 7;
        aRingTool.detectRingFromAtoms(allAtoms, allRings, dLev, nMaxRing);
        
        std::cout <<  "(2) Different tool. There are " 
                  << allRings.size() << " rings. They are: "
                  << std::endl;
        std::vector<RingDict> allRingsV;    
        for (std::map<std::string, std::vector<LIBMOL::RingDict> > ::iterator iR1=allRings.begin();
                    iR1 !=allRings.end(); iR1++)
        {
            std::cout << "(2)Ring representation " << iR1->first << std::endl;
            for (std::vector<RingDict>::iterator iR11=iR1->second.begin();
                        iR11 !=iR1->second.end(); iR11++)
            {
                allRingsV.push_back(*iR11);
                std::cout << "The ring consists of atoms: " << std::endl;
                for (std::vector<AtomDict>::iterator iAt1=iR11->atoms.begin();
                        iAt1 !=iR11->atoms.end(); iAt1++)
                {
                    std::cout << iAt1->id << std::endl;
                }
            }
                
            std::cout << std::endl;
            
        }
        
        if (allRingsV.size())
        {
            checkAndSetupPlanes(allRingsV, allPlanes, allAtoms);
        }
            
        aRingTool.setAtomsRingRepreS(allAtoms, allRingsV);
        
        for (int i=0; i < (int)allAtoms.size(); i++)
        {
            allAtoms[i].codClass = "";
            setAtomCodClassNameNew(allAtoms[i], allAtoms[i], dLev);
            std::cout <<std::endl << "For atom " << allAtoms[i].id << std::endl 
                      << "class is " << allAtoms[i].codClass << std::endl;            
        }
        
    }
    
    void CodClassify::codAtomClassifyNew2(int dLev)
    {
        ringTools aRingTool;
        allRings.clear();
        allPlanes.clear();
        int nMaxRing = 7;
        aRingTool.detectRingFromAtoms(allAtoms, allRings, dLev, nMaxRing);
        
        std::cout <<  "(2) Different tool. There are " 
                  << allRings.size() << " rings. They are: "
                  << std::endl;
        std::vector<RingDict> allRingsV;    
        for (std::map<std::string, std::vector<LIBMOL::RingDict> > ::iterator iR1=allRings.begin();
                    iR1 !=allRings.end(); iR1++)
        {
            std::cout << "(2)Ring representation " << iR1->first << std::endl;
            for (std::vector<RingDict>::iterator iR11=iR1->second.begin();
                        iR11 !=iR1->second.end(); iR11++)
            {
                allRingsV.push_back(*iR11);
                std::cout << "The ring consists of atoms: " << std::endl;
                for (std::vector<AtomDict>::iterator iAt1=iR11->atoms.begin();
                        iAt1 !=iR11->atoms.end(); iAt1++)
                {
                    std::cout << iAt1->id << std::endl;
                }
            }
                
            std::cout << std::endl;
            
        }
        
        if (allRingsV.size())
        {
            checkAndSetupPlanes(allRingsV, allPlanes, allAtoms);
            setAromaticBonds(allRingsV, allBonds);
        }
            
        aRingTool.setAtomsRingRepreS(allAtoms, allRingsV);
        
        for (int i=0; i < (int)allAtoms.size(); i++)
        {
            allAtoms[i].codClass = "";
            setAtomCodClassNameNew2(allAtoms[i], allAtoms[i], dLev);
            
            //std::cout <<std::endl << "For atom " << allAtoms[i].id << std::endl 
            //          << "class is " << allAtoms[i].codClass << std::endl;            
        }
        
        // merge into above late on
        for (std::vector<AtomDict>::iterator iAt=allAtoms.begin();
                iAt !=  allAtoms.end(); iAt++)
        {
           setSpecial3NBSymb2(iAt);
           //std::cout << "Now atom " << iAt->id << std::endl 
           //          << "class is " << iAt->codClass << std::endl; 
        }
        
        setAtomsNBSymb2();
        
        hashingAtoms2();
        
        
        
    }
    
    void CodClassify::codAtomClassifyNew3(int dLev)
    {
        ringTools aRingTool;
        allRings.clear();
        allPlanes.clear();
        int nMaxRing = 7;
        aRingTool.detectRingFromAtoms(allAtoms, allRings, dLev, nMaxRing);
        
        std::cout <<  "(2) Different tool. There are " 
                  << allRings.size() << " rings. They are: "
                  << std::endl;
        
        std::vector<RingDict> allRingsV;    
        for (std::map<std::string, std::vector<LIBMOL::RingDict> > ::iterator iR1=allRings.begin();
                    iR1 !=allRings.end(); iR1++)
        {
            std::cout << "(2)Ring representation " << iR1->first << std::endl;
            for (std::vector<RingDict>::iterator iR11=iR1->second.begin();
                        iR11 !=iR1->second.end(); iR11++)
            {
                allRingsV.push_back(*iR11);
                std::cout << "The ring consists of atoms: " << std::endl;
                for (std::vector<AtomDict>::iterator iAt1=iR11->atoms.begin();
                        iAt1 !=iR11->atoms.end(); iAt1++)
                {
                    std::cout << iAt1->id << std::endl;
                }
            }
                
            std::cout << std::endl;
            
        }
        
        if (allRingsV.size())
        {
            checkAndSetupPlanes(allRingsV, allPlanes, allAtoms);
        }
            
        aRingTool.setAtomsRingRepreS(allAtoms, allRingsV);
        
        for (int i=0; i < (int)allAtoms.size(); i++)
        {
            allAtoms[i].codClass = "";
            setAtomCodClassNameNew2(allAtoms[i], allAtoms[i], dLev);
            
            std::cout << std::endl << "For atom " << allAtoms[i].id << std::endl 
                      << "class is " << allAtoms[i].codClass << std::endl;            
        }
        
        for (std::vector<AtomDict>::iterator iAt=allAtoms.begin();
               iAt !=  allAtoms.end(); iAt++)
        {
            setSpecial3NBSymb(iAt);
            std::cout << "Now atom " << iAt->id << std::endl 
                      << "class is " << iAt->codClass << std::endl; 
        }
        
        
    }
    
    void CodClassify::getSmallFamily(std::string tInStr, NB1stFam& aNBFam)
    {
        std::vector<std::string> ch_list;
        ch_list.push_back(",");
        ch_list.push_back("x");

        std::string name_str = "";
        bool l_r = false;
        int  n_rep =1;
        
        for (int i =0; i < (int)tInStr.size(); i++)
        {
            if (std::isalpha(tInStr[i]))
            {
                char c;
                c=tInStr[i];
                if (std::toupper(c) ==tInStr[i])
                {
                    if (name_str.size() !=0)
                    {
                        if (aNBFam.name =="")
                        {
                            aNBFam.name = name_str;
                            n_rep =1;
                        }
                        else
                        {
                            for (int l=0; l < n_rep; l++)
                            {
                                aNBFam.NB2ndList.push_back(name_str);
                            }
                            n_rep =1;
                        }
                    }
                    name_str = tInStr[i];
                }
                else
                {
                    name_str= name_str +tInStr[i];
                }
            }
            else if (tInStr[i]=='[')
            {
                name_str= name_str +tInStr[i];
                l_r  = true;
            }
            else if (tInStr[i]==(']'))
            {
                name_str= name_str +tInStr[i];
                l_r  = false;
            }
            else if (std::find(ch_list.begin(), ch_list.end(), tInStr.substr(i,1)) !=ch_list.end())
            {
                name_str= name_str +tInStr[i];
            }
            else if (std::isdigit(tInStr[i]))
            {
                if (l_r)
                {
                    name_str= name_str +tInStr[i];
                }
                else
                {
                    n_rep = StrToInt(tInStr.substr(i,1));
                }
            }
        }
        
            
        if (aNBFam.name =="")
        {
            aNBFam.name = name_str;
        }
        else
        {
            for (int l=0;  l < n_rep; l++)
            {
                aNBFam.NB2ndList.push_back(name_str);
            }
        }
        
    }
    
    void CodClassify::codClassToAtom(ID  & tCodClass, AtomDict& tAt)
    {
        tCodClass=TrimSpaces(tCodClass);
        
        tAt.codClass = tCodClass;
        
        std::vector<std::string> atmStrs;
        
        StrTokenize(tCodClass, atmStrs, '(');
        
        int nRpt1 =0;
        
        if((int)atmStrs.size() !=0)
        {
            size_t f=atmStrs[0].find('[');
            if( f!=std::string::npos)
            {
                std::vector<ID> aPros;
                StrTokenize(atmStrs[0], aPros, '[');
                tAt.cChemType=TrimSpaces(aPros[0]);
                tAt.ringRep[aPros[1]] = 1;  // May need further actions 
            }
            else
            {
                tAt.cChemType=TrimSpaces(atmStrs[0]);
            }
            std::string aC = tAt.cChemType.substr(0,1);
            StrUpper(aC);
            tAt.chemType = tAt.cChemType.replace(0,1, aC);
            if(aC !=tAt.cChemType.substr(0,1))
            {
                tAt.bondingIdx = 2;
            }
            
            // Now two neighbor symbols
            
            for (int i=1; i <(int)atmStrs.size(); i++)
            {
                nRpt1 = 0;
                std::vector<ID> aNBgrp;
                StrTokenize(TrimSpaces(atmStrs[i]), aNBgrp, ')');
                if((int)aNBgrp.size()==2 and (int)aNBgrp[1].size() !=0)
                {     
                    nRpt1= StrToInt(aNBgrp[1]);
                }
                
                
                // get secondary neighbor atoms and the associated first NB atom
                int n2NB =0;
                
                
                ID  symb1 ="";
                ID  symb2 ="";
                ID  rep    ="";
                ID  achar  ="";        
                
                f=aNBgrp[0].find("]");
                if (f !=aNBgrp[0].npos)
                {
                    symb1   += aNBgrp[0].substr(0,(int)f+1);
                    symb2    = aNBgrp[0].substr((int)f+1);
                }
                else
                {
                    int iPos =0;
                    for(int j=0; j < (int)aNBgrp[0].size(); j++)
                    {
                        achar=aNBgrp[0].substr(j,1);
                        StrUpper(achar);
                        if(achar.compare(aNBgrp[0].substr(j,1))!=0)
                        {
                            symb1 +=aNBgrp[0].substr(j,1);
                            iPos+=1;
                        }
                        else
                        {
                            if ((int)symb1.size()==0)
                            {
                                symb1=aNBgrp[0].substr(j,1);
                                iPos+=1;
                            }
                            else
                            {
                                break;
                            }
                        }
                        
                    }
                    symb2 = aNBgrp[0].substr(iPos);
                }    
                    
                rep   = "";
                std::vector <ID> nb2list;
                for (int j=0; j < (int)symb2.size();  j++)
                {
                    ID achar = symb2.substr(j,1);
                    StrUpper(achar);
                    if (achar.compare(symb2.substr(j,1))!=0)
                    {
                        rep+=achar;
                    }
                    else if(isInt(achar))
                    {
                        int nchar = StrToInt(achar);
                        for(int k=0; k<nchar; k++)
                        {
                            nb2list.push_back(rep);
                        }
                            
                        rep = "";
                    }
                    else
                    {
                        if ((int)rep.size()==0)
                        {
                            rep =symb2[j];
                        }
                        else
                        {
                            nb2list.push_back(rep);
                            rep =symb2[j];
                        }
                    }
                }
                if((int)rep.size())
                {
                    nb2list.push_back(rep);
                }  
                
                n2NB = (int)nb2list.size()+1;
                symb1 +=("-"+IntToStr(n2NB)+":");
                if(nRpt1==0)
                {
                    tAt.codNBSymb +=symb1;
                    tAt.codNB2Symb +=(IntToStr(n2NB)+":");
                }
                else
                {
                    for (int i=0; i < nRpt1; i++)
                    {
                        tAt.codNBSymb  +=symb1;
                        tAt.codNB2Symb +=(IntToStr(n2NB)+":");
                    }
                }
            }
        }
    }
    
    
    void CodClassify::codClassToAtom2(ID  & tCodClass, AtomDict& tAt)
    {
        
        std::vector<NB1stFam>      allNBs;
        
        tCodClass=TrimSpaces(tCodClass);
        
        tAt.codClass = tCodClass;
        
        std::vector<std::string> tTwoP;
        std::string tMainSec;
        
        if (tCodClass.find('{') !=std::string::npos)
        {
            StrTokenize(tCodClass, tTwoP, '{');
            if (tTwoP.size()==2)
            {
                tMainSec = tTwoP[0];
                std::vector<std::string> ttTwoP;
                StrTokenize(tTwoP[1], ttTwoP, '}');
                tAt.codNB3Symb = ttTwoP[0];
                tAt.codAtmMain = tTwoP[0];
            }
            else
            {
                tAt.codAtmMain= tCodClass;
            }
        }
        else
        {
             tAt.codAtmMain = tCodClass;
        }
        //std::cout << "Atom class is " << tCodClass << std::endl;
        //std::cout << "Main section " << tAt.codAtmMain << std::endl;
        
        
        
        std::vector<std::string> atmStrs;
        
        StrTokenize(tAt.codAtmMain, atmStrs, '(');
        
        if((int)atmStrs.size() !=0)
        {
            tAt.codAtmRoot =  atmStrs[0];
            
            size_t f=atmStrs[0].find('[');
            if( f!=std::string::npos)
            {
                std::vector<ID> aPros;
                StrTokenize(atmStrs[0], aPros, '[');
                tAt.cChemType=TrimSpaces(aPros[0]);
                
                // tAt.ringRep[aPros[1]] = 1;  // May need further actions 
            }
            else
            {
                tAt.cChemType=TrimSpaces(atmStrs[0]);
            }
            tAt.chemType = tAt.cChemType;
            
            
            
            // Now Set all level of atom hierarchical symbols
            
            if (atmStrs.size() < 2)
            {
                std::cout << tAt.id << " has atom class symbol: "
                          << tAt.codClass << ", which has no NB, wrong string ? " 
                          << std::endl;
                exit(1);
            }
            
            
            //std::cout << "Atom itself " << tAt.codAtmRoot << std::endl;
 
            //std::cout << "Split Neighbor symbols: "  << std::endl;
            
            for (unsigned i=1; i < atmStrs.size(); i++)
            {
                
                std::string tS = TrimSpaces(atmStrs[i]);

                NB1stFam aNBFam;
                std::vector<std::string> NB1;
                std::vector<std::string> NB1Main;
                StrTokenize(tS, NB1, ')' );
                if (NB1.size() > 1)
                {
                    aNBFam.repN=  StrToInt(NB1[1]);
                    if (aNBFam.repN==0)
                    {
                        aNBFam.repN = 1;
                    }
                }
                else
                {
                    aNBFam.repN = 1;
                }


                std::string tS1=TrimSpaces(NB1[0]);
                //std::cout << "One NB string :  " << tS1 << std::endl;
                //std::cout << "Unit          :  " << aNBFam.repN << std::endl;

                getSmallFamily(tS1, aNBFam);
                /*
                for (int j=0; j < aNBFam.repN; j++)
                {
                    std::cout << "1st NB is " << aNBFam.name << std::endl;
    
                    std::cout << "This 1st NB has " << aNBFam.NB2ndList.size() << " second NB " << std::endl;
                    
                    if (aNBFam.NB2ndList.size())
                    {
                        //std::cout << "Those 2NB are: " << std::endl;
    
                        for (unsigned k=0; k < aNBFam.NB2ndList.size(); k++)
                        {
                            std::cout << aNBFam.NB2ndList[k] << std::endl;
                        }
                        
                        
                    }
                }
                */
                
                allNBs.push_back(aNBFam);
            }
        }
        
        if ( tAt.codNB3Symb.size() !=0)
        {
            std::cout << "The 3rd NB part is " << tAt.codNB3Symb << std::endl;
        }
        
        
        for (std::vector<NB1stFam>::iterator iNB=allNBs.begin();
                iNB != allNBs.end(); iNB++)
        {
            for (int j=0; j < iNB->repN; j++)
            {
                ID sN=IntToStr((int)iNB->NB2ndList.size()+1);
                tAt.codNBSymb  += (iNB->name + "-" + sN + ":");
                tAt.codNB2Symb +=(sN+":");
            }
        }
        
        std::cout << "Summary: " << std::endl
                  << "Atom class " << tAt.codClass << std::endl
                  << "Atom class main section " << tAt.codAtmMain << std::endl
                  << "Level 2 atom NB symbol "  << tAt.codNB2Symb << std::endl
                  << "Level 3 atom NB symbol "  << tAt.codNBSymb  << std::endl;
        if (tAt.codNB3Symb.size())
        {
            std::cout << "The 3rd NB section is "   << tAt.codNB3Symb << std::endl;
        }
                  
    }
    
    
    void CodClassify::codClassToAtomAng(ID  & tCodClass, AtomDict& tAt)
    {
        tCodClass=TrimSpaces(tCodClass);
        
        // tAt.codClass  = tCodClass;
        tAt.codNBSymb = "";
        tAt.codNB2Symb ="";
        
        std::vector<std::string> atmStrs;
        
        StrTokenize(tCodClass, atmStrs, '(');
        
        int nRpt1 =0;
        
        if((int)atmStrs.size() !=0)
        {
            size_t f=atmStrs[0].find('[');
            if( f!=std::string::npos)
            {
                std::vector<ID> aPros;
                StrTokenize(atmStrs[0], aPros, '[');
                tAt.cChemType=TrimSpaces(aPros[0]);
                tAt.ringRep[aPros[1]] = 1;  // May need further actions 
            }
            else
            {
                tAt.cChemType=TrimSpaces(atmStrs[0]);
            }
            std::string aC = tAt.cChemType.substr(0,1);
            StrUpper(aC);
            tAt.chemType = tAt.cChemType.replace(0,1, aC);
            if(aC !=tAt.cChemType.substr(0,1))
            {
                tAt.bondingIdx = 2;
            }
            
            // Now two neighbor symbols
            
            
            for (int i=1; i <(int)atmStrs.size(); i++)
            {
                nRpt1 = 0;
                std::vector<ID> aNBgrp;
                StrTokenize(TrimSpaces(atmStrs[i]), aNBgrp, ')');
                if((int)aNBgrp.size()==2 and (int)aNBgrp[1].size() !=0)
                {     
                    nRpt1= StrToInt(aNBgrp[1]);
                }
                
                
                // get secondary neighbor atoms and the associated first NB atom
                int n2NB =0;
                
                
                ID  symb1 ="";
                ID  symb2 ="";
                ID  rep    ="";
                ID  achar  ="";        
                
                f=aNBgrp[0].find("]");
                if (f !=aNBgrp[0].npos)
                {
                    symb1   += aNBgrp[0].substr(0,(int)f+1);
                    symb2    = aNBgrp[0].substr((int)f+1);
                }
                else
                {
                    int iPos =0;
                    for(int j=0; j < (int)aNBgrp[0].size(); j++)
                    {
                        achar=aNBgrp[0].substr(j,1);
                        StrUpper(achar);
                        if(achar.compare(aNBgrp[0].substr(j,1))!=0)
                        {
                            symb1 +=aNBgrp[0].substr(j,1);
                            iPos+=1;
                        }
                        else
                        {
                            if ((int)symb1.size()==0)
                            {
                                symb1=aNBgrp[0].substr(j,1);
                                iPos+=1;
                            }
                            else
                            {
                                break;
                            }
                        }
                        
                    }
                    symb2 = aNBgrp[0].substr(iPos);
                }    
                    
                rep   = "";
                std::vector <ID> nb2list;
                for (int j=0; j < (int)symb2.size();  j++)
                {
                    ID achar = symb2.substr(j,1);
                    StrUpper(achar);
                    if (achar.compare(symb2.substr(j,1))!=0)
                    {
                        rep+=achar;
                    }
                    else if(isInt(achar))
                    {
                        int nchar = StrToInt(achar);
                        for(int k=0; k<nchar; k++)
                        {
                            nb2list.push_back(rep);
                        }
                            
                        rep = "";
                    }
                    else
                    {
                        if ((int)rep.size()==0)
                        {
                            rep =symb2[j];
                        }
                        else
                        {
                            nb2list.push_back(rep);
                            rep =symb2[j];
                        }
                    }
                }
                if((int)rep.size())
                {
                    nb2list.push_back(rep);
                }  
                
                ID tNBAtmChemType;
                
                std::vector<std::string> tNBAtmStrs;
                StrTokenize(symb1, tNBAtmStrs, '[');
                
                tNBAtmChemType= tNBAtmStrs[0];
                
                std::vector<ID> tMetab;
                initMetalTab(tMetab);
                
                if (!isMetal(tMetab,tNBAtmChemType))
                {
                    n2NB = (int)nb2list.size()+1;
                    symb1 +=("-"+IntToStr(n2NB)+":");
                    if(nRpt1==0)
                    {
                        tAt.codNBSymb +=symb1;
                        tAt.codNB2Symb +=(IntToStr(n2NB)+":");
                    }
                    else
                    {
                        for (int i=0; i < nRpt1; i++)
                        {
                            tAt.codNBSymb +=symb1;
                            tAt.codNB2Symb +=(IntToStr(n2NB)+":");
                        }
                    }
                }
            }
        }
    }
    
    void CodClassify::codNBProps(std::vector<ID> tarStrs, std::vector<ID>& tCTs,
                                 std::vector<int> & tNB2s, 
                                 std::vector<int> & tRis, std::vector<int>& tPls)
    {
        
        for (int i=0; i < (int)tarStrs.size(); i++)
        {
            std::vector<ID>  tS;
            StrTokenize(tarStrs[i], tS, '-');
            tNB2s.push_back(StrToInt(tS[1]));
            size_t f=tS[0].find('[');
            if( f!=std::string::npos)
            {   
                
                tRis.push_back(1);
                std::vector<ID> aPros;
                StrTokenize(tS[0], aPros, '[');
                ID c1 = aPros[0].substr(0,1);
                StrUpper(c1);
                if (c1.compare(aPros[0].substr(0,1)) !=0)
                {
                    tPls.push_back(1);
                }
                else
                {
                    tPls.push_back(0);
                }
                if (aPros[0].size()==0)
                {
                    tCTs.push_back(c1);
                }
                else
                {
                    tCTs.push_back(c1+aPros[0].substr(1));
                }
            }
            else
            {
                tRis.push_back(0);
                ID c1 = tS[0].substr(0,1);
                StrUpper(c1);
                if (c1.compare(tS[0].substr(0,1)) !=0)
                {
                    tPls.push_back(1);
                }
                else
                {
                    tPls.push_back(0);
                }
                if (tS[0].size()==0)
                {
                    tCTs.push_back(c1);
                }
                else
                {
                    tCTs.push_back(c1+tS[0].substr(1));
                }
            }
        }
    }
  
    int CodClassify::codAtomsDist(std::vector<ID> tarStrs, std::vector<ID> tNBs, 
                                  int tLev)
    {
        // In atom hashing the following features are listed as required matches:
        // (1) chemical types of root atoms, 
        // (2) if they are in planes 
        // (3) if they are in rings 
        // (4) numbers of the first neighbor atoms
        // 
        // The distances contributed by different features of atoms
        // (1) differences in numbers of the secondary neighbor atom 
        // (2) differences in chemical types of the first NB atoms 
         
        int allDist = 0;
        std::vector<ID> tarCTs,bCTs;
        std::vector<int> tarNB2s, tarRis, tarPls, bNB2s, bRis, bPls;

        codNBProps(tarStrs, tarCTs, tarNB2s, tarRis, tarPls);
        codNBProps(tNBs, bCTs, bNB2s, bRis, bPls);
        if (tLev ==1)
        {
            allDist=codAtomsChemTypeDist(tarCTs, tarRis, tarPls,
                                         bCTs, bRis, bPls);
        }
        else if (tLev==2)
        {
            allDist = codAtomsChemTypeDist(tarCTs, tarRis, tarPls,
                                           bCTs, bRis, bPls)
                    + codAtomsNumNBDist(tarNB2s, bNB2s);
        }
        
        return allDist;
    }
    
    int CodClassify::codAtomsChemTypeDist(std::vector<ID>  & tarCTs,
                                          std::vector<int>  & tarRis, 
                                          std::vector<int>  & tarPls,
                                          std::vector<ID>  & bCTs,
                                          std::vector<int>  & bRis,
                                          std::vector<int>  & bPls)
    {
        int tDist =0;
            
        int tUnit1 = 20;  // for row
        int tUnit2 = 100; // for group
        int tUnit3 = 50;  // for plane and ring
            
        for (int i=0; i <(int)tarCTs.size(); i++)
        {
            int i11 = pPeriodictable->elements[tarCTs[i]]["row"];
            int i12 = pPeriodictable->elements[tarCTs[i]]["group"];
            int i21 = pPeriodictable->elements[bCTs[i]]["row"];
            int i22 = pPeriodictable->elements[bCTs[i]]["group"];
            tDist+=(std::abs(i11-i21)*tUnit1);
            tDist+=(std::abs(i12-i22)*tUnit2);
            tDist+=(std::abs(tarRis[i]-bRis[i])*tUnit3);
            tDist+=(std::abs(tarPls[i]-bPls[i])*tUnit3);
        }
        
        return tDist;
        
    }
    
    int CodClassify::codAtomsNumNBDist(std::vector<int> & tNB2s,
                                       std::vector<int> & bNB2s)
    {
        int tDist = 0;
        int tUnit = 150;
        
        for (int i=0; i <(int)tNB2s.size(); i++)
        {
            tDist+=(std::abs(tNB2s[i]-bNB2s[i])*tUnit);
        }
        
        return tDist;
    }
    
    // Ring related 
    void CodClassify::ringDetecting()
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
            /*
                std::cout << "-----------------" << std::endl;
                std::cout << "starting atom  " << iA->id << std::endl;
                std::cout << "-----------------" << std::endl;
            */
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
            
        
        }
       
    }
    
    void CodClassify::ringDetecting2()
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
            /*
                std::cout << "-----------------" << std::endl;
                std::cout << "starting atom  " << iA->id << std::endl;
                std::cout << "-----------------" << std::endl;
            */
            //if((int)tempIDs.size() !=0)
            //{
            //    tempIDs.clear();
            //}
   
            if (!iA->isMetal )
            {
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
                checkOnePathSec2(*iA, maxSize, iA, preSeriNum,  startLev, atomsSeen, atomsInPath);
            
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
            }
        }
       
       
    }
    
    
    void CodClassify::checkOnePathSec(std::vector<AtomDict> & seenAtoms,
                                      std::map<int, ID>     & seenIDs,
                                      AtomDict              & curAto,
                                      int                     iMax,
                                      std::vector<AtomDict>::iterator iAto)

    {
        int tSizeSeen   = (int)seenAtoms.size();
        int tPreSeriNum;
        std::string tPreId;
        
        if (tSizeSeen >=2)
        {
            tPreSeriNum= seenAtoms[tSizeSeen-2].seriNum;
            tPreId     = seenAtoms[tSizeSeen-2].id;
        }
        else
        {
            tPreSeriNum = 0;
            tPreId      = "N/a";
        }
        if ( tSizeSeen <iMax )
        {
            // loop bonded atoms in the current atom
            for (std::vector<int>::iterator tNBA=curAto.connAtoms.begin();
                    tNBA != curAto.connAtoms.end(); tNBA++)
            {
                int tSeriNum = allAtoms[*tNBA].seriNum;
                
                
                // should not go back to the atom searched in the last step
                if (tSeriNum != tPreSeriNum)
                {
                        
                    //std::cout << "Current seen atom size: " << (int)seenAtoms.size() << std::endl;
                    //std::cout << "beginning  atom: " << iAto->id << std::endl;
                    //std::cout << "Previous atom: " << tPreId << std::endl;
                    //std::cout << "Current  atom: " << curAto.id << std::endl;
                    //std::cout << "Pass Next atom " << allAtoms[*tNBA].id << std::endl;
                    //std::cout << tPreSeriNum << " : " << curAto.seriNum << " : " 
                    //              << tSeriNum << std::endl;
                      
                    
                    std::map<int, ID>::iterator tLI = seenIDs.find(tSeriNum);
                    if (tLI!=seenIDs.end() && (int)seenIDs.size() !=1)
                    {
                        // find the atom in the seen-atom list 
                        // Check if this atom is a Nachbarpunkte atom
                        if (tSeriNum != iAto->seriNum)
                        {
                            // is a Nachbarpunkte atom, stop this path
                            
                            //   std::cout << " is a Nachbarpunkte atom, stop this path " << std::endl;
                            break;
                        }
                        else
                        {
                            // find a ring 
                           
                            
                            //   std::cout << "Find a ring. It contains  " << std::endl;
                            
                            // sort the atoms in the seenAtom vector and 
                            // check if this ring is already found (the same 
                            // ring of the same atom should be found twice at
                            // least because of the walking algorithm. 
                            std::list<std::string> tAllIds;
                            
                            for (int iSee=1; iSee < (int)seenAtoms.size(); iSee++)
                            {
                                tAllIds.push_back(seenAtoms[iSee].id);
                                // tRing.atoms.push_back(*tAt);
                                //    std::cout << seenAtoms[iSee].id << std::endl;
                                
                            }
                            tAllIds.sort(compareNoCase);
                            std::string tRepStr(seenAtoms[0].id);
                            for (std::list<std::string>::iterator iAI =tAllIds.begin();
                                    iAI != tAllIds.end(); iAI++)
                            {
                                tRepStr.append(*iAI);
                            }
                            
                           
                            //    std::cout << "The ring representation (after sorting) is : "
                            //               << tRepStr << std::endl;
                           
                            if (seenAtoms[0].ringRep.find(tRepStr) 
                                    == seenAtoms[0].ringRep.end())
                            {
                                seenAtoms[0].ringRep.insert(std::pair<std::string,int>(tRepStr,1));
                                
                                //   std::cout << "This ring has been added to the ring list of the atom"
                                //              << std::endl;
                                
                                Ring tRing;
                                for (std::vector<AtomDict>::iterator iSAtom =seenAtoms.begin();
                                        iSAtom != seenAtoms.end(); iSAtom++)
                                {
                                    tRing.atoms.push_back(*iSAtom);
                                }
                                
                                
                                // std::cout << "Begin atom is " << iAto->id << std::endl;
                                iAto->ringRep[tRepStr] = (int)tRing.atoms.size();
                            }
                        }
                    }
                    else 
                    {
                        // This is a new atom and append the atom to seen-atom list
                        // and call function checkOnePathSec() recursively
                        seenAtoms.push_back(allAtoms[*tNBA]);
                        seenIDs.insert(std::pair<int, ID>(allAtoms[*tNBA].seriNum,
                                                          allAtoms[*tNBA].id));
                        
                        checkOnePathSec(seenAtoms, seenIDs, allAtoms[*tNBA], iMax, iAto);
                        
                        // if last call of checkOnePathSec does not find a ring.
                        // remove the current atom from the seen-list and loop 
                        seenAtoms.pop_back();
                        seenIDs.erase(allAtoms[*tNBA].seriNum);
                    }
                }
            }
        }
    }
        
    
    // Version 2
    
    void CodClassify::checkOnePathSec(AtomDict                & curAto,
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
    
    void CodClassify::checkOnePathSec2(AtomDict                & curAto,
                                      int                        iMax,
                                      std::vector<AtomDict>::iterator iOriAto,
                                      int                        SeriNumPreAto,  
                                      int                        curLev,
                                      std::map<int, ID>       &  seenAtomIDs,
                                      std::map<int, ID>       &  atomIDsInPath)
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
                        std::list<std::string> tAllSeris;
                        std::vector <AtomDict> ttAtoms;
                        for (std::map<int, ID>::iterator iSee = atomIDsInPath.begin();
                                iSee != atomIDsInPath.end(); iSee++)
                        {
                            //if (iSee->first != iOriAto->seriNum)
                            //{
                            tAllSeris.push_back( IntToStr(iSee->first));
                            tAllIds.push_back(iSee->second);
                            
                            //int posIdx = atomPosition(iSee->second);
                            ttAtoms.push_back(allAtoms[iSee->first]);
                            
                            // tRing.atoms.push_back(*tAt);
                            //std::cout << iSee->second << std::endl;   
                            //}
                        }
                        
                        RingDict aRingDict(ttAtoms);                 
                        tAllSeris.sort(compareNoCase);
                        tAllIds.sort(compareNoCase);
                        //std::string tRepStr(iOriAto->id);
                        std::string tRepStr;
                        std::string tRepSeri;
                        for (std::list<std::string>::iterator iAI =tAllIds.begin();
                                    iAI != tAllIds.end(); iAI++)
                        {
                            tRepStr.append(*iAI);
                        }
                        
                        int nRS =0;
                        for (std::list<std::string>::iterator iAS =tAllSeris.begin();
                                    iAS != tAllSeris.end(); iAS++)
                        {
                            if (nRS==0)
                            {
                                tRepSeri.append(*iAS);
                            }
                            else
                            {
                                tRepSeri.append("_" + *iAS);
                            }
                            nRS++;
                        }
                            
                        iOriAto->ringRep[tRepStr] = (int)atomIDsInPath.size();
                        
                        aRingDict.rep  = tRepStr;
                        aRingDict.sRep = tRepSeri;
                        
                        std::map<ID, std::vector<RingDict> >::iterator iFindRing=allRings.find(tRepSeri);
                        if (iFindRing == allRings.end())
                        {
                            for (std::map<int, ID>::iterator iSee = atomIDsInPath.begin();
                                iSee != atomIDsInPath.end(); iSee++)
                            {
                                // int posIdx = atomPosition(iSee->second);
                                allAtoms[iSee->first].inRings.push_back((int)allRings.size());
                            }
                            allRings[tRepSeri].push_back(aRingDict);
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
                        if (SeriNumPreAto != allAtoms[*tNBA].seriNum && ! allAtoms[*tNBA].isMetal)
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
                            checkOnePathSec2(allAtoms[*tNBA], iMax, iOriAto, tPreSeriNum, 
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
    
    void CodClassify::setAtomCodClassName(AtomDict & tAtom,
                                          AtomDict & tOriAtom,
                                          int tLev)
    {
        
        if (tLev==1)
        {
            tAtom.codClass = "";
            tAtom.codClass.append(tAtom.chemType);
            outRingSec(tAtom);
            
            std::string tStr;
            std::list<std::string> tStrList, tStrList1;
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
            
            sortMap  tCMap;
            std::vector<sortMap> tCVec;
            
            for (std::map<ID, int>::iterator tM=comps.begin();
                   tM !=comps.end(); tM++)
            {
                tCMap.key = tM->first;
                tCMap.val = tM->second;
                tCVec.push_back(tCMap);
            }
            
            std::sort(tCVec.begin(),tCVec.end(), desSortMapKey);
            
            for (std::vector<sortMap>::iterator iMa=tCVec.begin();
                    iMa !=tCVec.end(); iMa++)
            {
                std::string s1, s2;
                s1 = iMa->key + IntToStr(iMa->val);
                for (int i=0; i < iMa->val; i++)
                {
                    s2.append(iMa->key);
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
            
            /*
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
            */
            
            // tStrList.sort(compareNoCase2);
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
                AtomDict aNBAtom(allAtoms[*tNBA]);
                setAtomCodClassName(aNBAtom, tOriAtom, lowLev);
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
                
                if(tIdMap.find(aNBAtom.codClass) !=tIdMap.end())
                {
                    tIdMap[aNBAtom.codClass]++;
                    
                }
                else
                {
                    tIdMap[aNBAtom.codClass] = 1;
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
                if (iV->val ==1)
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
    
    void CodClassify::setAtomCodClassName2(AtomDict & tAtom,
                                           AtomDict & tOriAtom,
                                           int tLev)
    {
        
        if (tLev==1)
        {
            tAtom.codClass = "";
            tAtom.codClass.append(tAtom.chemType);
            outRingSec2(tAtom);
            
            std::string tStr;
            std::list<std::string> tStrList, tStrList1;
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
            
            sortMap  tCMap;
            std::vector<sortMap> tCVec;
            
            for (std::map<ID, int>::iterator tM=comps.begin();
                   tM !=comps.end(); tM++)
            {
                tCMap.key = tM->first;
                tCMap.val = tM->second;
                tCVec.push_back(tCMap);
            }
            
            std::sort(tCVec.begin(),tCVec.end(), desSortMapKey);
            
            for (std::vector<sortMap>::iterator iMa=tCVec.begin();
                    iMa !=tCVec.end(); iMa++)
            {
                std::string s1, s2;
                s1 = iMa->key + IntToStr(iMa->val);
                for (int i=0; i < iMa->val; i++)
                {
                    s2.append(iMa->key);
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
            
            /*
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
            */
            
            // tStrList.sort(compareNoCase2);
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
            std::map<std::string, std::vector<int> > tIdMap;
            for (std::vector<int>::iterator tNBA=tAtom.connAtoms.begin();
                    tNBA != tAtom.connAtoms.end(); tNBA++)
            {
                AtomDict aNBAtom(allAtoms[*tNBA]);
                setAtomCodClassName(aNBAtom, tOriAtom, lowLev);
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
                
                if(tIdMap.find(aNBAtom.codClass) !=tIdMap.end())
                {
                    tIdMap[aNBAtom.codClass][0]++;
                    
                }
                else
                {
                    tIdMap[aNBAtom.codClass].push_back(1);
                    tIdMap[aNBAtom.codClass].push_back((int)aNBAtom.connAtoms.size());
                }
                 
            }
            
            //sortMap  tSMap;
            //std::vector<sortMap> tVec;
            
            
            std::vector<sortMap2> tVec;
            
            for (std::map<std::string, std::vector<int> >::iterator tM=tIdMap.begin();
                   tM !=tIdMap.end(); tM++)
            {
                struct sortMap2  tSMap2;
                tSMap2.key = tM->first;
                tSMap2.val = tM->second[0];
                tSMap2.nNB = tM->second[1];
                tVec.push_back(tSMap2);
            }
            
            std::sort(tVec.begin(),tVec.end(), desSortMapKey2);
            
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
            for(std::vector<sortMap2>::iterator iV=tVec.begin();
                    iV !=tVec.end(); iV++)
            {
                if (iV->val ==1)
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
    
    void CodClassify::setAtomCodClassNameNew(AtomDict & tAtom,
                                           AtomDict & tOriAtom,
                                           int tLev)
    {
        
        if (tLev==1)
        {
            tAtom.codClass = "";
            tAtom.codClass.append(tAtom.chemType);
            outRingSecNew(tAtom);
            
            std::string tStr;
            std::list<std::string> tStrList, tStrList1;
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
            
            sortMap  tCMap;
            std::vector<sortMap> tCVec;
            
            for (std::map<ID, int>::iterator tM=comps.begin();
                   tM !=comps.end(); tM++)
            {
                tCMap.key = tM->first;
                tCMap.val = tM->second;
                tCVec.push_back(tCMap);
            }
            
            std::sort(tCVec.begin(),tCVec.end(), desSortMapKey);
            
            for (std::vector<sortMap>::iterator iMa=tCVec.begin();
                    iMa !=tCVec.end(); iMa++)
            {
                std::string s1, s2;
                s1 = iMa->key + IntToStr(iMa->val);
                for (int i=0; i < iMa->val; i++)
                {
                    s2.append(iMa->key);
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
            
            /*
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
            */
            
            // tStrList.sort(compareNoCase2);
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
            outRingSecNew(tAtom);
            //std::cout << "Atom " << tAtom.id << " its COD ring section " 
            //          <<  tAtom.codClass << std::endl;
            
            int lowLev = tLev - 1;
            std::map<std::string, std::vector<int> > tIdMap;
            for (std::vector<int>::iterator tNBA=tAtom.connAtoms.begin();
                    tNBA != tAtom.connAtoms.end(); tNBA++)
            {
                AtomDict aNBAtom(allAtoms[*tNBA]);
                setAtomCodClassNameNew(aNBAtom, tOriAtom, lowLev);
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
                
                if(tIdMap.find(aNBAtom.codClass) !=tIdMap.end())
                {
                    tIdMap[aNBAtom.codClass][0]++;
                    
                }
                else
                {
                    tIdMap[aNBAtom.codClass].push_back(1);
                    tIdMap[aNBAtom.codClass].push_back((int)aNBAtom.connAtoms.size());
                }
                 
            }
            
            //sortMap  tSMap;
            //std::vector<sortMap> tVec;
            
            
            std::vector<sortMap2> tVec;
            
            for (std::map<std::string, std::vector<int> >::iterator tM=tIdMap.begin();
                   tM !=tIdMap.end(); tM++)
            {
                struct sortMap2  tSMap2;
                tSMap2.key = tM->first;
                tSMap2.val = tM->second[0];
                tSMap2.nNB = tM->second[1];
                tVec.push_back(tSMap2);
            }
            
            std::sort(tVec.begin(),tVec.end(), desSortMapKey2);
            
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
            for(std::vector<sortMap2>::iterator iV=tVec.begin();
                    iV !=tVec.end(); iV++)
            {
                if (iV->val ==1)
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
    
    void CodClassify::setAtomCodClassNameNew2(AtomDict & tAtom,
                                             AtomDict & tOriAtom,
                                             int tLev)
    {
        
        if (tLev==1)
        {
            tAtom.codClass = "";
            tAtom.codClass.append(tAtom.chemType);
            outRingSecNew(tAtom);
            
            std::string tStr;
            std::list<std::string> tStrList, tStrList1;
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
                    std::string t2NBType;
                    t2NBType.append(allAtoms[*tNBAtom].chemType);
                    outRingSecNew2(t2NBType, allAtoms[*tNBAtom]);
                    if(comps.find(t2NBType) != comps.end())
                    {
                        comps[t2NBType] += 1;
                    }
                    else
                    {
                        comps[t2NBType] = 1; 
                    }
                }
            }
            
            sortMap  tCMap;
            std::vector<sortMap> tCVec;
            
            for (std::map<ID, int>::iterator tM=comps.begin();
                   tM !=comps.end(); tM++)
            {
                tCMap.key = tM->first;
                tCMap.val = tM->second;
                tCVec.push_back(tCMap);
            }
            
            std::sort(tCVec.begin(),tCVec.end(), desSortMapKey);
            
            for (std::vector<sortMap>::iterator iMa=tCVec.begin();
                    iMa !=tCVec.end(); iMa++)
            {
                std::string s1, s2;
                s1 = iMa->key + IntToStr(iMa->val);
                for (int i=0; i < iMa->val; i++)
                {
                    s2.append(iMa->key);
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
            
            /*
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
            */
            
            // tStrList.sort(compareNoCase2);
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
            outRingSecNew(tAtom);
            //std::cout << "Atom " << tAtom.id << " its COD ring section " 
            //          <<  tAtom.codClass << std::endl;
            
            int lowLev = tLev - 1;
            std::map<std::string, std::vector<int> > tIdMap;
            for (std::vector<int>::iterator tNBA=tAtom.connAtoms.begin();
                    tNBA != tAtom.connAtoms.end(); tNBA++)
            {
                AtomDict aNBAtom(allAtoms[*tNBA]);
                setAtomCodClassNameNew2(aNBAtom, tOriAtom, lowLev);
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
                
                if(tIdMap.find(aNBAtom.codClass) !=tIdMap.end())
                {
                    tIdMap[aNBAtom.codClass][0]++;
                    
                }
                else
                {
                    tIdMap[aNBAtom.codClass].push_back(1);
                    tIdMap[aNBAtom.codClass].push_back((int)aNBAtom.connAtoms.size());
                }
                 
            }
            
            //sortMap  tSMap;
            //std::vector<sortMap> tVec;
            
            
            std::vector<sortMap2> tVec;
            
            for (std::map<std::string, std::vector<int> >::iterator tM=tIdMap.begin();
                   tM !=tIdMap.end(); tM++)
            {
                struct sortMap2  tSMap2;
                tSMap2.key = tM->first;
                tSMap2.val = tM->second[0];
                tSMap2.nNB = tM->second[1];
                tVec.push_back(tSMap2);
            }
            
            std::sort(tVec.begin(),tVec.end(), desSortMapKey2);
            
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
            for(std::vector<sortMap2>::iterator iV=tVec.begin();
                    iV !=tVec.end(); iV++)
            {
                if (iV->val ==1)
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
    
    
    void CodClassify::outRingSec(AtomDict &tAtom)
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
    
    void CodClassify::outRingSec2(AtomDict &tAtom)
    {
        int numRings = (int)tAtom.ringRepBySeriNum.size();
       
        if (numRings)
        {
            
            std::map<int, int> sizeMap;
            
            
            for (std::map<std::string, int>::iterator iMR=tAtom.ringRepBySeriNum.begin();
                    iMR != tAtom.ringRepBySeriNum.end(); iMR++)
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
    
    
    void CodClassify::outRingSecNew(AtomDict &tAtom)
    {
        int numRings = (int)tAtom.ringRepS.size();
        //std::cout << "tAtom.ringRepS.size() " << tAtom.ringRepS.size() << std::endl;
        
        if (numRings)
        {
            
            std::map<std::string, int> sizeMap;
            
            
            for (std::map<std::string, std::string>::iterator iMR=tAtom.ringRepS.begin();
                    iMR != tAtom.ringRepS.end(); iMR++)
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
            for (std::map<ID, int>::iterator iSMa=sizeMap.begin();
                    iSMa != sizeMap.end(); iSMa++)
            {
                // std::string tSize = IntToStr(iSMa->first);
                std::string tSize(iSMa->first);
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
            
            //std::cout << "ring section atom " << tAtom.id << tAtom.codClass << std::endl;
        }
    }
    
    void CodClassify::outRingSecNew2(std::string & tAtmCodStr,
                                     AtomDict &tAtom)
    {
        int numRings = (int)tAtom.ringRepS.size();
        //std::cout << "tAtom.ringRepS.size() " << tAtom.ringRepS.size() << std::endl;
        
        if (numRings)
        {
            
            std::map<std::string, int> sizeMap;
            
            
            for (std::map<std::string, std::string>::iterator iMR=tAtom.ringRepS.begin();
                    iMR != tAtom.ringRepS.end(); iMR++)
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
            
            tAtmCodStr.append("[");
            int i =0;
            int j = (int)sizeMap.size();
            for (std::map<ID, int>::iterator iSMa=sizeMap.begin();
                    iSMa != sizeMap.end(); iSMa++)
            {
                // std::string tSize = IntToStr(iSMa->first);
                std::string tSize(iSMa->first);
                std::string tNum  = IntToStr(iSMa->second);
                
                if(iSMa->second >= 3)
                {
                    tAtmCodStr.append(tNum + "x" + tSize);
                }
                else if (iSMa->second==2)
                {
                    tAtmCodStr.append( tSize + "," + tSize);
                }   
                else if (iSMa->second==1)
                {
                    tAtmCodStr.append(tSize);
                }
                       
                
                if(i != j-1)
                {
                    tAtmCodStr.append(",");
                }
                else
                {
                    tAtmCodStr.append("]");
                }
                    
                i++;
            }
            
            //std::cout << "ring section atom " << tAtom.id << tAtom.codClass << std::endl;
        }
    }
    
    
    std::string CodClassify::outRingSecStr(AtomDict &tAtom)
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
    
    void CodClassify::detectPlaneGroups()
    {
        
        groupOrgAtomsToPlanes();
        groupMetAndLigandAtomsToPlanes();
   
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
    
    void CodClassify::setAtomsBondingAndChiralCenter()
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
            //std::cout << "Atom " << iAt->id << std::endl
            //        <<  " connect to  " << t_len << std::endl;
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
                else if (iAt->connAtoms.size()==1)
                {
                    iAt->bondingIdx = 1;
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
            // std::cout << "its chiralIdx " << iAt->chiralIdx << std::endl;
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
                    if (iAt->parCharge ==0.0)
                    {
                        bool l_sp2 = false;
                        for (std::vector<int>::iterator iCA=iAt->connAtoms.begin();
                                 iCA != iAt->connAtoms.end(); iCA++)
                        {
                            if(allAtoms[*iCA].bondingIdx == 2)
                            {
                                l_sp2 = true;
                                break;
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
                    else if (iAt->parCharge ==1.0)
                    {
                        iAt->chiralIdx  =  0;
                        iAt->bondingIdx =  2;
                    }
                } 
            }
            //std::cout << "Again atom " << iAt->id << " its chiralIdx " 
            //          << iAt->chiralIdx << std::endl;
        }
        
        // Further check if a chiral center is a real one
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA != allAtoms.end(); iA++)
        {
            if (iA->chiralIdx !=0)
            {
                std::vector<ID> chirRAtms;
                for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                        iNB != iA->connAtoms.end(); iNB++)
                {
                    std::size_t tFind = allAtoms[*iNB].chemType.find("H");
                    if (tFind !=std::string::npos)
                    {
                        chirRAtms.push_back(allAtoms[*iNB].id);
                    }
                }
                if ((int)chirRAtms.size() >1 && (int)iA->connAtoms.size() <=4)
                {
                    iA->chiralIdx = 0;
                }
            }
            
            
        }
        /*
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
                //int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3;
                }
                else if (t_len ==3)
                {
                    iAt->chiralIdx  = -1;
                    iAt->bondingIdx = 2;
                } 
            }
            else if (iAt->chemType.compare("N")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3;  
                }
                //else if (t_len==3) // temp 
                //{   // should do on the next round when all NB atoms are set
                //    iAt->chiralIdx  = -1;
                //    iAt->bondingIdx =  2;
                // }
                else if (t_len ==2)
                {
                    iAt->chiralIdx  = -1;
                    iAt->bondingIdx =  2;
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
            else if (iAt->chemType.compare("SI")==0 
                    || iAt->chemType.compare("P")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3; 
                }
                else if (t_len==3)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 2; 
                }
            }
            else if (iAt->chemType.compare("S")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                if(t_len==4 || t_len==3)
                {
                    iAt->chiralIdx  = 1;
                    iAt->bondingIdx = 3; 
                }
            }
        }
        
        // more conditions 
        
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {            
            if (iAt->chemType.compare("N")==0 || iAt->chemType.compare("B")==0)
            {
                // int t_len = (int)iAt->connAtoms.size();
                int t_len =0;
                for (std::vector<int>::iterator iConn=iAt->connAtoms.begin();
                     iConn !=iAt->connAtoms.end(); iConn++)
                {
                    if(!allAtoms[*iConn].isMetal)
                    {
                        t_len++;
                    }
                }
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
                        iAt->chiralIdx  = -1;
                        iAt->bondingIdx =  2;
                    }
                    else
                    {
                        iAt->chiralIdx  =  1;
                        iAt->bondingIdx =  3;
                    }
                } 
            }
        }
        
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA != allAtoms.end(); iA++)
        {
            if (iA->chiralIdx ==1)
            {
                std::vector<ID> chirRAtms;
                for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                        iNB != iA->connAtoms.end(); iNB++)
                {
                    std::size_t tFind = allAtoms[*iNB].chemType.find("H");
                    if (tFind !=std::string::npos)
                    {
                        chirRAtms.push_back(allAtoms[*iNB].id);
                    }
                }
                if ((int)chirRAtms.size() >1 && (int)iA->connAtoms.size() <=4)
                {
                    iA->chiralIdx = 0;
                }
            }
        }*/
        
        
        // No need for the third round, those could be defined in 
    
        // Check
        /*
        std::cout << "Chiral and plane feather for atoms in the system" 
                  << std::endl;
        
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {
            if (iAt->chiralIdx == -1)
            {
                std::cout << "Atom " << iAt->id << " may be  in a negative chiral center " 
                        << std::endl;
            }
            else if (iAt->chiralIdx == 1)
            {
                std::cout << "Atom " << iAt->id << " may be in a positive chiral center "
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
        */
        
       
    }
    
    void CodClassify::groupOrgAtomsToPlanes()
    {
        std::vector<PlaneDict> smalPls;
        
        // Find the smallest planes
        setSmallestPLs(smalPls);
        
        // Merge the above planes to the large plane groups
        mergeLargePLGroups(smalPls);
        
    }
    

    void CodClassify::setSmallestPLs(std::vector<PlaneDict>& tSmaPls)
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
        /*
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
        */
        
    }
    
    void CodClassify::mergeLargePLGroups(std::vector<PlaneDict>& tSmaPls)
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

    void CodClassify::groupMetAndLigandAtomsToPlanes()
    {
        for (int i=0; i < (int)allAtoms.size(); i++)
        {
            if (allAtoms[i].isMetal)
            {
                if (allAtoms[i].metalGeo.compare("BENT")==0 
                    || allAtoms[i].metalGeo.compare("LINEAR") ==0
                    || allAtoms[i].metalGeo.compare("T-SHAPE")==0
                    || allAtoms[i].metalGeo.compare("TRIGONAL-PLANAR")==0
                    || allAtoms[i].metalGeo.compare("SQUARE-PLANAR")==0)
                {
                    PlaneDict aPl;
                    aPl.archID = allAtoms[i].id;
                    aPl.archPos = i;
                    aPl.atoms[allAtoms[i].id] = i;
                    for (std::vector<int>::iterator iN=allAtoms[i].connAtoms.begin();
                            iN !=allAtoms[i].connAtoms.end(); iN++)
                    {
                        aPl.atoms[allAtoms[*iN].id] = *iN;
                    }
                    allPlanes.push_back(aPl);
                }
            }
        }
    }
    
    bool CodClassify::isInSameRing(PlaneDict & tP1, PlaneDict & tP2)
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
    
    bool CodClassify::furtherM(std::vector<int> &tV1, std::vector<int> &tV2)
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
    
    void CodClassify::readTablesForHashing(std::map<int, ID> & tDigitKeys,
                                           std::map<int, int> & tLinkedHA)
    {
        
        for (int i =0; i < wSize+1; i++)
        {
            tDigitKeys[i] = NullString;
            tLinkedHA[i]  = -1;
        }
        
        //std::string clibMonDir(std::getenv("CLIBD_MON"));
        //std::string inFileLHAName = clibMonDir + "allOrgLinkedHashCode.table";
        //std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
        //std::string inFileLHAName = clibMonDir + "/tables/allOrgLinkedHashCode.table";
        
        
        std::string inFileLHAName = libmolTabDir  + "/allOrgLinkedHashCode.table";    
        std::ifstream inFileLHA;
        inFileLHA.open(inFileLHAName.c_str(), std::ios::in);
        //inFileLHA.open("/Users/flong/COD/New_EXP/Current/derivedData/allOrgLinkedHashCode.table",
        //               std::ios::in);
        
        if(inFileLHA.is_open())
        {
            std::string tRecord="";
            while(!inFileLHA.eof())
            {
                std::getline(inFileLHA, tRecord);
                tRecord = TrimSpaces(tRecord);
                std::vector<std::string> tBuf;
                StrTokenize(tRecord, tBuf);
                
                if((int)tBuf.size() > 2)
                {
                    int aHA = StrToInt(tBuf[0]);
                    tDigitKeys[aHA] = tBuf[1];
                    tLinkedHA[aHA]  = StrToInt(tBuf[2]);
                }
            }
            
            inFileLHA.close();
        }
        else
        {
            std::cout << "can not find the tables at " << inFileLHAName << std::endl;
            exit(1);
        }
    }
    
    void CodClassify::readTablesForHashing2(std::map<int, ID> & tDigitKeys,
                                            std::map<int, int> & tLinkedHA)
    {
        
        for (int i =0; i < wSize+1; i++)
        {
            tDigitKeys[i] = NullString;
            tLinkedHA[i]  = -1;
        }
        
        //std::string clibMonDir(std::getenv("CLIBD_MON"));
        //std::string inFileLHAName = clibMonDir + "allOrgLinkedHashCode.table";
        //std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
        //std::string inFileLHAName = clibMonDir + "/tables/allOrgLinkedHashCode.table";
        
        
        std::string inFileLHAName = libmolTabDir  + "/allOrgLinkedHashCode.table";  
        // std::string inFileLHAName = "/Applications/ccp4-6.5/share/acedrg/tables_DB3/allOrgLinkedHashCode.table";
        std::ifstream inFileLHA;
        inFileLHA.open(inFileLHAName.c_str(), std::ios::in);
        //inFileLHA.open("/Users/flong/COD/New_EXP/Current/derivedData/allOrgLinkedHashCode.table",
        //               std::ios::in);
        
        if(inFileLHA.is_open())
        {
            std::string tRecord="";
            while(!inFileLHA.eof())
            {
                std::getline(inFileLHA, tRecord);
                tRecord = TrimSpaces(tRecord);
                std::vector<std::string> tBuf;
                StrTokenize(tRecord, tBuf);
                
                if((int)tBuf.size() > 2)
                {
                    int aHA = StrToInt(tBuf[0]);
                    tDigitKeys[aHA] = tBuf[1];
                    tLinkedHA[aHA]  = StrToInt(tBuf[2]);
                }
            }
            
            inFileLHA.close();
        }
        else
        {
            std::cout << "can not find the tables at " << inFileLHAName << std::endl;
            exit(1);
        }
    }
    
    void CodClassify::hashingAtoms()
    {
        
        std::map<int, ID>    aDigitKeys;
        std::map<int, int>   aLinkedHA;
        
        readTablesForHashing(aDigitKeys, aLinkedHA);
        
        std::vector<int> aPrimTab;
        initPrimeTab(aPrimTab, libmolTabDir);
        // std::cout << aPrimTab.size() << std::endl;
     
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {
            std::map<ID, std::map<std::string, int> >::iterator iFind;
            iFind = pPeriodictable->elements.find(iAt->chemType);
            if (iFind !=pPeriodictable->elements.end())
            {
                int iMin = iAt->getMinRing();
                //std::cout << "Atom " << iAt->codClass << " has min ring "
                //        << iMin << std::endl;
               
                int d2, d3, d4, d5;
                /*
                if (iAt->bondingIdx == 2)
                {
                    d1 = 1;
                }
                else
                {
                    d1 = 0;
                }
                
                // if ((int)iAt->ringRep.size())
                {
                    d2 = 3;
                }
                else
                {
                    d2 = 2; 
                }
                */
                if (iMin == 0)
                {
                    d2 = 1;
                }
                else if (iMin ==3)
                {
                    d2 = 2;
                }
                else if (iMin==4)
                {
                    d2 = 3;
                }
                else if (iMin==5) 
                {
                    d2 = 4;
                }
                else if (iMin==6)
                {
                    d2 = 5;
                }
                else
                {
                    d2 = 6;
                }
                    
                
                
                d3 = 7 + (int)iAt->connAtoms.size();
                
                d4 = 15 + pPeriodictable->elements[iAt->chemType]["row"];
                
                d5 = 23 + pPeriodictable->elements[iAt->chemType]["group"];
                
                /*
                std::cout << "chemType " << iAt->chemType << std::endl
                          << "codClass " << iAt->codClass << std::endl;
                std::cout << "row " << pPeriodictable->elements[iAt->chemType]["row"] << std::endl;
                std::cout << "group "<< pPeriodictable->elements[iAt->chemType]["group"] << std::endl;
                
                std::cout << " d1 " << d1 << " d2 " << d2 << " d3 " << d3 
                         << " d4 " << d4  << " d5 " << d5 << std::endl;
                
                */
                
                int aPrim = aPrimTab[d2]*aPrimTab[d3]*aPrimTab[d4]*aPrimTab[d5];               
                
                int psedoHA   = aPrim%wSize;
                ID  footPrint = IntToStr(d2) + "_" + IntToStr(d3) 
                                + "_" + IntToStr(d4) + "_" + IntToStr(d5);
 
                
                std::map<int, ID>::iterator  dkFinder=aDigitKeys.find(psedoHA);
                
                if (dkFinder != aDigitKeys.end())
                {
                    if (aDigitKeys[psedoHA].compare(footPrint)==0)
                    {
                        iAt->hashingValue = psedoHA;
                    }
                    else
                    {
                        bool lCont = true;
                        while(lCont)
                        {
                            std::map<int, int>::iterator lkFinder=aLinkedHA.find(psedoHA);
                            if (lkFinder !=aLinkedHA.end())
                            {
                                int aHA = aLinkedHA[psedoHA];
                                if (aHA==-1)
                                {
                                    std::cout << "No hash code in DB for atom " 
                                            << iAt->codClass << std::endl;
                                    lCont = false;
                                }
                                else
                                {
                                    if(aDigitKeys[aHA]==footPrint)
                                    {
                                        iAt->hashingValue = aHA;
                                        lCont = false;
                                    }
                                    else
                                    {
                                        psedoHA=aHA;
                                    }
                                }  
                            }
                            else 
                            {
                                std::cout << "No hash code in DB for atom " 
                                          << iAt->codClass << std::endl;
                                lCont = false;
                            }
                        }
                    }
                }
                else
                {
                    std::cout << "No footprint in DB for atom " 
                              << iAt->codClass << std::endl;
                }
            }
            else
            {
                std::cout << "No element ID in DB for atom " 
                          << iAt->codClass << std::endl;
            }
            
           //std::cout << "Atom " << iAt->id << " has hashing code : " 
           //          << iAt->hashingValue << std::endl;
        }   
        
    }
    
    void CodClassify::hashingAtoms2()
    {
        
        std::map<int, ID>    aDigitKeys;
        std::map<int, int>   aLinkedHA;
        
        readTablesForHashing2(aDigitKeys, aLinkedHA);
        
        std::vector<int> aPrimTab;
        // libmolTabDir ="/Applications/ccp4-6.5/share/acedrg/tables_DB3";
        initPrimeTab(aPrimTab, libmolTabDir);
        
        
        // std::cout << aPrimTab.size() << std::endl;
     
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {
            std::map<ID, std::map<std::string, int> >::iterator iFind;
            iFind = pPeriodictable->elements.find(iAt->chemType);
            if (iFind !=pPeriodictable->elements.end())
            {
                int iMin = iAt->getMinRing2();
               
                //std::cout << "Atom " << iAt->codClass << " has min ring "
                //          << iMin << std::endl;
               
                int d1,d2, d3, d4, d5;
                /*
                if (iAt->bondingIdx == 2)
                {
                    d1 = 1;
                }
                else
                {
                    d1 = 0;
                }
                
                // if ((int)iAt->ringRep.size())
                {
                    d2 = 3;
                }
                else
                {
                    d2 = 2; 
                }
                */
                
                if (iAt->baseRingProp["aroma"].compare("y") ==0)
                {
                    d1 = 1;
                }
                else
                {
                    d1 = 0;
                }
                
                if (iMin == 0)
                {
                    d2 = 2;
                }
                else if (iMin ==3)
                {
                    d2 = 3;
                }
                else if (iMin==4)
                {
                    d2 = 4;
                }
                else if (iMin==5) 
                {
                    d2 = 5;
                }
                else if (iMin==6)
                {
                    d2 = 6;
                }
                else
                {
                    d2 = 7;
                }
                    
                
                
                d3 = 8 + (int)iAt->connAtoms.size();
                
                d4 = 16 + pPeriodictable->elements[iAt->chemType]["row"];
                
                d5 = 24 + pPeriodictable->elements[iAt->chemType]["group"];
                
                
                std::cout << "chemType " << iAt->chemType << std::endl
                          << "codClass " << iAt->codClass << std::endl;
                std::cout << "row " << pPeriodictable->elements[iAt->chemType]["row"] << std::endl;
                std::cout << "group "<< pPeriodictable->elements[iAt->chemType]["group"] << std::endl;
                
                std::cout << " d1 " << d1 << " d2 " << d2 << " d3 " << d3 
                          << " d4 " << d4  << " d5 " << d5 << std::endl;
                
                
                
                int aPrim = aPrimTab[d1]*aPrimTab[d2]*aPrimTab[d3]*aPrimTab[d4]*aPrimTab[d5];               
                
                int psedoHA   = aPrim%wSize;
                ID  footPrint = IntToStr(d1) + "_" + IntToStr(d2) + "_" + IntToStr(d3) 
                                + "_" + IntToStr(d4) + "_" + IntToStr(d5);
 
                //std::cout << "psedoHA " << psedoHA << std::endl;
                //std::cout << "footPrint " << footPrint << std::endl;
                
                std::map<int, ID>::iterator  dkFinder=aDigitKeys.find(psedoHA);
                
                if (dkFinder != aDigitKeys.end())
                {   
                    if (aDigitKeys[psedoHA].compare(footPrint)==0)
                    {
                        iAt->hashingValue = psedoHA;
                    }
                    else
                    {
                        bool lCont = true;
                        while(lCont)
                        {
                            std::map<int, int>::iterator lkFinder=aLinkedHA.find(psedoHA);
                            if (lkFinder !=aLinkedHA.end())
                            {
                                int aHA = aLinkedHA[psedoHA];
                                if (aHA==-1)
                                {
                                    std::cout << "No hash code in DB for atom " 
                                            << iAt->codClass << std::endl;
                                    lCont = false;
                                }
                                else
                                {
                                    if(aDigitKeys[aHA]==footPrint)
                                    {
                                        iAt->hashingValue = aHA;
                                        lCont = false;
                                    }
                                    else
                                    {
                                        psedoHA=aHA;
                                    }
                                }  
                            }
                            else 
                            {
                                std::cout << "No hash code in DB for atom " 
                                          << iAt->codClass << std::endl;
                                lCont = false;
                            }
                        }
                    }
                }
                else
                {
                    std::cout << "No footprint in DB for atom " 
                              << iAt->codClass << std::endl;
                }
            }
            else
            {
                std::cout << "No element ID in DB for atom " 
                          << iAt->codClass << std::endl;
            }
            
            std::cout << "Atom " << iAt->id << " has hashing code : " 
                      << iAt->hashingValue << std::endl;
        }   
        
    }
    
    
    
    void CodClassify::setAtomsNBSymb()
    {
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {
            
            AtomDict tAtm;
            codClassToAtom(iAt->codClass, tAtm);
            iAt->codNBSymb = tAtm.codNBSymb;
            iAt->codNB2Symb = tAtm.codNB2Symb;
            /*
            std::list<ID> tList1;
            
            for (std::vector<int>::iterator iNBidx=iAt->connAtoms.begin();
                    iNBidx !=iAt->connAtoms.end(); iNBidx++)
            {
                std::vector<std::string> tGrpStrs;
                StrTokenize(allAtoms[*iNBidx].codClass, tGrpStrs, '(');
                if ((int)tGrpStrs.size()!=0)
                {
                    std::string numNB = IntToStr((int)allAtoms[*iNBidx].connAtoms.size())+ ":";
                    tList1.push_back(tGrpStrs[0]+ "-" + numNB);
                }
            }
            
            tList1.sort(compareNoCase);
            for (std::list<std::string>::iterator iAI =tList1.begin();
                 iAI != tList1.end(); iAI++)
            {
                iAt->codNBSymb+=(*iAI);
                std::vector<ID> tS1,tS2;
                StrTokenize(*iAI, tS1, '-');
                StrTokenize(tS1[1], tS2, ':');
                iAt->codNB2Symb+=(TrimSpaces(tS2[0])+":");
            }
            */
            
            //std::cout << iAt->id << std::endl;
            //std::cout << iAt->codClass   << std::endl;
            //std::cout << iAt->codNBSymb  << std::endl;
            //std::cout << iAt->codNB2Symb << std::endl;
        }
          
       
        
    }
    
    
    void CodClassify::setAtomsNBSymb2()
    {
        for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                iAt != allAtoms.end(); iAt++)
        {
            
            AtomDict tAtm;
            // std::cout << "The atom serial number " << iAt->seriNum << std::endl; 
            codClassToAtom2(iAt->codClass, tAtm);
            iAt->codAtmMain = tAtm.codAtmMain;
            iAt->codNBSymb  = tAtm.codNBSymb;
            iAt->codNB2Symb = tAtm.codNB2Symb;
            iAt->codNB3Symb = tAtm.codNB3Symb;
            
            /*
            std::list<ID> tList1;
            
            for (std::vector<int>::iterator iNBidx=iAt->connAtoms.begin();
                    iNBidx !=iAt->connAtoms.end(); iNBidx++)
            {
                std::vector<std::string> tGrpStrs;
                StrTokenize(allAtoms[*iNBidx].codClass, tGrpStrs, '(');
                if ((int)tGrpStrs.size()!=0)
                {
                    std::string numNB = IntToStr((int)allAtoms[*iNBidx].connAtoms.size())+ ":";
                    tList1.push_back(tGrpStrs[0]+ "-" + numNB);
                }
            }
            
            tList1.sort(compareNoCase);
            for (std::list<std::string>::iterator iAI =tList1.begin();
                 iAI != tList1.end(); iAI++)
            {
                iAt->codNBSymb+=(*iAI);
                std::vector<ID> tS1,tS2;
                StrTokenize(*iAI, tS1, '-');
                StrTokenize(tS1[1], tS2, ':');
                iAt->codNB2Symb+=(TrimSpaces(tS2[0])+":");
            }
            */
            
            //std::cout << iAt->id << std::endl;
            //std::cout << iAt->codClass   << std::endl;
            //std::cout << iAt->codNBSymb  << std::endl;
            //std::cout << iAt->codNB2Symb << std::endl;
        }
          
       
        
    }
    
    
    
    // Initiate  a set of dummy atoms for the atom tree 
    // and torsion angle calculations
    
    void CodClassify::initDummyAtoms()
    {

    }
    
    // Construct a MST for the system
    
    void CodClassify::setAtomsMST()
    {
        std::vector<int> atmsInT, atmsOutT;
        // supplement variables
        std::map<int, std::map<int, REAL> >  allDs;
        for (int i=0; i < (int)allAtoms.size(); i++)
        {
            for (std::vector<int>::iterator iNB=allAtoms[i].connAtoms.begin();
                    iNB != allAtoms[i].connAtoms.end(); iNB++)
            {
                allDs[i][*iNB] = distanceV(allAtoms[i].coords, allAtoms[*iNB].coords);
            }
            if ((int)atmsInT.size() ==0)
            {
                // The starting atom of the tree should not be H and not in rings
                if (allAtoms[i].id != "H" && (int)allAtoms[i].ringRep.size() ==0)
                {
                    // start point of the tree
                    atmsInT.push_back(i);
                    allAtoms[i].tree["parent"].push_back(-1);
                    std::cout << "Atom " << allAtoms[i].id << " is the starting atom " << std::endl;
                }
                else
                {
                    atmsOutT.push_back(i);
                }
            }
            else
            {
                atmsOutT.push_back(i);
            }
        }
        
        
        while ((int)atmsOutT.size() !=0)
        {
            std::cout << (int)atmsOutT.size() << " atoms out of tree yet " << std::endl;
            int iP=-1, iC=-1; 
            REAL dMin =2000.0;
            std::vector<int>::iterator iFind;
            for (std::vector<int>::iterator iTA=atmsInT.begin();
                    iTA !=atmsInT.end(); iTA++)
            {
                for (std::map<int, REAL>::iterator iNB=allDs[*iTA].begin();
                        iNB !=allDs[*iTA].end(); iNB++)
                {
                    std::vector<int>::iterator tFind = 
                       std::find(atmsOutT.begin(), atmsOutT.end(), iNB->first);
                    if (iNB->second < dMin && 
                         tFind!=atmsOutT.end())
                    {
                        iP=*iTA;
                        iC=iNB->first;
                        dMin = iNB->second;
                        iFind = tFind;
                    }
                }
            }
            
            // std::cout << "Parent " << iP << " child " << iC << std::endl;
            if (iP != -1 && iC !=-1)
            {
                // get one vertex of the tree
                allAtoms[iP].tree["children"].push_back(iC);
                allAtoms[iC].tree["parent"].push_back(iP);
                atmsInT.push_back(iC);
                std::cout << "Atom " << allAtoms[iP].id << " has a child " << allAtoms[iC].id 
                        << std::endl;
           
                std::cout << "Atom " << allAtoms[iC].id << " has parent " << allAtoms[iP].id << std::endl;
                // std::cout << " number of parent " << (int)allAtoms[iP].tree["parent"].size() << std::endl;
                
                if ((int)allAtoms[iC].tree["parent"].size() > 1)
                {
                    std::cout << "Atom " << allAtoms[iC].id << " has "
                            << (int)allAtoms[iC].tree["parent"].size() << std::endl
                            << "They are \t";
                    for (std::vector<int>::iterator iPA=allAtoms[iC].tree["parent"].begin();
                            iPA != allAtoms[iC].tree["parent"].begin(); iPA++)
                    {
                        std::cout << *iPA << "\t";
                    }
                    std::cout << std::endl;
                    
                    exit(1);            
                }
                /*
                else
                {
                    std::cout << "Atom " << allAtoms[iP].id << " has "
                            << (int)allAtoms[iP].tree["children"].size() 
                            << " children. They are \n";
                    for (std::vector<int>::iterator iPA=allAtoms[iP].tree["children"].begin();
                            iPA != allAtoms[iP].tree["children"].end(); iPA++)
                    {
                        std::cout << *iPA << "\t";
                    }
                    std::cout << std::endl;
                    
                }
                 */
                       
               
                // remove it for out of tree list and allDs
               
                atmsOutT.erase(iFind);
                
                allDs[iP].erase(iC);
                
            }
            else
            {
                std::cout << "why no further vertex found " << std::endl;
                std::cout << "Current out-of-tree atoms are " << std::endl;
                for (std::vector<int>::iterator iNTA=atmsOutT.begin();
                        iNTA !=atmsOutT.end(); iNTA++)
                {
                    std::cout << *iNTA << "\t";
                }
                std::cout << std::endl;          
            }
        }
       
        
        // Output tree for a look
        std::cout <<"******* Tree structure for all atoms " << std::endl;
        
        for (std::vector<AtomDict>::iterator iAt=allAtoms.begin(); 
                iAt != allAtoms.end(); iAt++)
        {
            
            
            if (iAt->tree.find("children") !=iAt->tree.end())
            {
                
                ID tPar;
                if (iAt->tree["parent"][0] !=-1)
                {
                    tPar = allAtoms[iAt->tree["parent"][0]].id; 
                }
                else
                {
                    tPar = "BEGIN";
                }
                    
                
                for (std::vector<int>::iterator iNA=iAt->tree["children"].begin();
                         iNA != iAt->tree["children"].end(); iNA++)
                {
                    std::cout << iAt->id << "\t"  << tPar << "\t"
                          << allAtoms[*iNA].id << std::endl;
                }
                
            }
            else
            {
                std::cout << iAt->id << "\t" << allAtoms[iAt->tree["parent"][0]].id 
                            << "\t" << "END" << std::endl; 
            }
        }
        
        
    }
    
    // Bond related 
    /*
    void CodClassify::setOrgBondHeadHashList()
    {
        //std::string clibMonDir(std::getenv("CLIBD_MON"));
        // std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
        for (std::vector<BondDict>::iterator iB = allBonds.begin();
                iB != allBonds.end(); iB++)
        {
            
            //std::string fRoot = clibMonDir + "allOrgBondTables/";
            std::string fRoot = libmolTabDir  + "allOrgBondTables/";
            int iSmall = 2*wSize;
            for (std::map<ID, int>::iterator iA=iB->fullAtoms.begin();
                    iA != iB->fullAtoms.end(); iA++)
            {
                if (allAtoms[iA->second].hashingValue < iSmall)
                {
                    iSmall = allAtoms[iA->second].hashingValue;
                }
            }
            std::map<int, ID>::iterator iFind = codOrgBondFiles.find(iSmall);
            if (iFind ==codOrgBondFiles.end())
            {
                codOrgBondFiles[iSmall] = fRoot + IntToStr(iSmall) + ".table";
            }
            
        }
        
        std::cout << "Should read following files :" << std::endl;
        for (std::map<int, ID>::iterator iBF = codOrgBondFiles.begin();
                iBF !=codOrgBondFiles.end(); iBF++)
        {
            std::cout << iBF->second << std::endl;
        } 
                
    }
    */
    
    void CodClassify::setOrgBondHeadHashList()
    {
        
        std::map<ID, std::map<ID, ID> >  allBoIdx;
        
        //std::string clibMonDir(std::getenv("CLIBD_MON"));
        //std::string fRoot = clibMonDir + "allOrgBondTables/";
        // std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
        
        std::string fRoot = libmolTabDir  + "/allOrgBondTables/";
        std::string fIdx  = fRoot + "bond_idx.table";
        std::ifstream codBondIdxFile(fIdx.c_str());
        if (codBondIdxFile.is_open())
        {
            std::string tRecord="";
            while(!codBondIdxFile.eof())
            {
                std::getline(codBondIdxFile, tRecord);    
                tRecord = TrimSpaces(tRecord);
                std::vector<std::string> tBuf;
                StrTokenize(tRecord, tBuf);
                if ((int)tBuf.size() ==3)
                {
                    allBoIdx[tBuf[0]][tBuf[1]] = tBuf[2];
                }
            }
            codBondIdxFile.close();
        }
        
                
        for (std::vector<BondDict>::iterator iBo = allBonds.begin();
                iBo != allBonds.end(); iBo++)
        {
            //std::string fRoot = "/Users/flong/COD/New_EXP/Current/derivedData/allOrgBondTables/";
            ID ha0, ha1;
            
            if (allAtoms[iBo->atomsIdx[0]].hashingValue <= allAtoms[iBo->atomsIdx[1]].hashingValue)
            {
                ha0 = IntToStr(allAtoms[iBo->atomsIdx[0]].hashingValue);
                ha1 = IntToStr(allAtoms[iBo->atomsIdx[1]].hashingValue);
            }
            else
            {
                ha0 = IntToStr(allAtoms[iBo->atomsIdx[1]].hashingValue);
                ha1 = IntToStr(allAtoms[iBo->atomsIdx[0]].hashingValue);
            }
                
            ID haNum;
            
            if (allBoIdx.find(ha0) != allBoIdx.end())
            {
                if (allBoIdx[ha0].find(ha1) != allBoIdx[ha0].end())
                {
                    haNum = allBoIdx[ha0][ha1];
                    std::map<ID, ID>::iterator iFind = codOrgBondFiles2.find(haNum);
                    if (iFind ==codOrgBondFiles2.end())
                    {
                        codOrgBondFiles2[haNum] = fRoot + haNum + ".table";
                    }
                }
            }
        }
       
        /*
        std::cout << "Bonds are in the following files :" << std::endl;
        for (std::map<ID, ID>::iterator iBF = codOrgBondFiles2.begin();
                iBF !=codOrgBondFiles2.end(); iBF++)
        {
            std::cout << iBF->second << std::endl;
        } 
        */
       
    }
    
    void CodClassify::setOrgBondHeadHashList2()
    {
        
        std::map<ID, std::map<ID, ID> >  allBoIdx;
        
        //std::string clibMonDir(std::getenv("CLIBD_MON"));
        //std::string fRoot = clibMonDir + "allOrgBondTables/";
        // std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
        //libmolTabDir ="/Applications/ccp4-6.5/share/acedrg/tables_DB3";
        std::string fRoot = libmolTabDir  + "/allOrgBondTables/";
        std::string fIdx  = fRoot + "bond_idx.table";
        std::ifstream codBondIdxFile(fIdx.c_str());
        
        if (codBondIdxFile.is_open())
        {
            std::string tRecord="";
            while(!codBondIdxFile.eof())
            {
                std::getline(codBondIdxFile, tRecord);    
                tRecord = TrimSpaces(tRecord);
                std::vector<std::string> tBuf;
                StrTokenize(tRecord, tBuf);
                if ((int)tBuf.size() ==3)
                {
                    allBoIdx[tBuf[0]][tBuf[1]] = tBuf[2];
                    // std::cout << "allBoIdx[" << tBuf[0] << "][" << tBuf[1] << "] =" << tBuf[2] << std::endl;
                }
            }
            codBondIdxFile.close();
        }
        
                
        for (std::vector<BondDict>::iterator iBo = allBonds.begin();
                iBo != allBonds.end(); iBo++)
        {
            
            ID ha0, ha1;
            //std::string fRoot = "/Users/flong/COD/New_EXP/Current/derivedData/allOrgBondTables/";
            /*
           
            if (tHa0 >= tHa1)
            {
                ha0 = tHa0;
                ha1 = tHa1;
            }
            else
            {
                ha0 =tHa1;
                ha1 =tHa0;
             }
             */
            
            if (allAtoms[iBo->atomsIdx[0]].hashingValue <= allAtoms[iBo->atomsIdx[1]].hashingValue)
            {
                ha0 = IntToStr(allAtoms[iBo->atomsIdx[0]].hashingValue);
                ha1 = IntToStr(allAtoms[iBo->atomsIdx[1]].hashingValue);
            }
            else
            {
                ha0 = IntToStr(allAtoms[iBo->atomsIdx[1]].hashingValue);
                ha1 = IntToStr(allAtoms[iBo->atomsIdx[0]].hashingValue);
            }
              
            ID haNum;
            
            if (allBoIdx.find(ha0) != allBoIdx.end())
            {
                if (allBoIdx[ha0].find(ha1) != allBoIdx[ha0].end())
                {
                    haNum = allBoIdx[ha0][ha1];
                    std::map<ID, ID>::iterator iFind = codOrgBondFiles2.find(haNum);
                    if (iFind ==codOrgBondFiles2.end())
                    {
                        codOrgBondFiles2[haNum] = fRoot + haNum + ".table";
                    }
                }
            }
        }
        
       
        
        std::cout << "Bonds are in the following files :" << std::endl;
        for (std::map<ID, ID>::iterator iBF = codOrgBondFiles2.begin();
                iBF !=codOrgBondFiles2.end(); iBF++)
        {
            std::cout << iBF->second << std::endl;
        } 
        
        
    }
 
    /*
    void  CodClassify::groupCodOrgBonds()
    {
        
        setOrgBondHeadHashList();
        time_t tStart, tEnd;
        std::time (&tStart);
        std::cout << "Clustering COD org bonds started at " << std::ctime(&tStart);
        int nline = 0;
        for (std::map<int, ID>::iterator iBF=codOrgBondFiles.begin();
                    iBF !=codOrgBondFiles.end(); iBF++)
        {
            try
            {
                // should be something like std::string tNewCodBondFileName(clibMonDir + "/list/bonds.txt");
                //std::string clibMonDir(std::getenv("CLIBD_MON"));
                //std::string codBondFileName = clibMonDir + "/allOrgBonds.table";
                //std::ifstream codBondFile(codBondFileName.c_str());
            
                std::ifstream codBondFile(iBF->second.c_str()); 
                if(codBondFile.is_open())
                {
                    std::string tRecord="";
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
                            // std::cout << "ha1 " << ha1 << " ha2 " << ha2 << std::endl;
                        
                        
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
                            aBond.value      = StrToReal(tBuf[8]);
                            aBond.sigValue   = StrToReal(tBuf[9]);
                            if(aBond.sigValue < 0.01)
                            {
                                aBond.sigValue = 0.01;
                            }
                            aBond.numCodValues = StrToInt(tBuf[10]);
                            //std::cout << aBond.atomsCodClasses[0] << std::endl;
                            //std::cout << aBond.atomsCodClasses[1] << std::endl;
                            //std::cout << nline << std::endl;
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
        }
      
        std::time(&tEnd);
        std::cout << "Clustering COD org bonds finished at " << std::ctime(&tEnd);
        REAL tDiff;
        tDiff = std::difftime(tEnd,tStart);
        std::cout  << "it takes " << std::setprecision(3) <<tDiff 
                   << " seconds to finish group COD bonds " << std::endl;
    }
    
    */
    
    
    void  CodClassify::groupCodOrgBonds()
    {
        setOrgBondHeadHashList();
        time_t tStart, tEnd;
       
        std::time (&tStart);
        std::cout << "Clustering COD org bonds started at " << std::ctime(&tStart);
        
        int nline = 0;
        for (std::map<ID, ID>::iterator iBF=codOrgBondFiles2.begin();
                    iBF !=codOrgBondFiles2.end(); iBF++)
        {
            try
            {
                // should be something like std::string tNewCodBondFileName(clibMonDir + "/list/bonds.txt");
                //std::string clibMonDir(std::getenv("CLIBD_MON"));
                //std::string codBondFileName = clibMonDir + "/allOrgBonds.table";
                //std::ifstream codBondFile(codBondFileName.c_str());
            
                std::ifstream codBondFile(iBF->second.c_str()); 
                if(codBondFile.is_open())
                {
                    std::string tRecord="";
                    while(!codBondFile.eof())
                    {
                        std::getline(codBondFile, tRecord);
                        
                        tRecord = TrimSpaces(tRecord);
                        std::vector<std::string> tBuf;
                        StrTokenize(tRecord, tBuf);
                        
                        if ((int)tBuf.size() ==20)
                        {
                            int ha1, ha2;
                        
                            ha1 = StrToInt(tBuf[0]);
                            ha2 = StrToInt(tBuf[1]);
                            // std::cout << "ha1 " << ha1 << " ha2 " << ha2 << std::endl;
                        
                            /*
                            for (int i=2; i < (int)tBuf.size(); i++)
                            {
                                tBuf[i] = TrimSpaces(tBuf[i]);
                            }
                            */
                            
                            allDictBondsIdx[ha1][ha2][tBuf[2]][tBuf[3]][tBuf[4]][tBuf[5]][tBuf[6]][tBuf[7]]= nline;  
                            
                            BondDict aBond;
                            aBond.seriNum = nline;
                            aBond.atomsHashingCodes.push_back(ha1);
                            aBond.atomsHashingCodes.push_back(ha2);
                            aBond.atomsNB2Rep.push_back(tBuf[2]);
                            aBond.atomsNB2Rep.push_back(tBuf[3]);
                            aBond.atomsNBRep.push_back(tBuf[4]);
                            aBond.atomsNBRep.push_back(tBuf[5]);
                            aBond.atomsCodClasses.push_back(tBuf[6]);
                            aBond.atomsCodClasses.push_back(tBuf[7]);
                            aBond.value      = StrToReal(tBuf[8]);
                            aBond.sigValue   = StrToReal(tBuf[9]);
                            //aBond.valueP     = StrToReal(tBuf[11]);
                            //aBond.sigValueP  = StrToReal(tBuf[12]);
                            if(aBond.sigValue < 0.01)
                            {
                                aBond.sigValue = 0.01;
                            }
                            aBond.numCodValues  = StrToInt(tBuf[10]);
                            //aBond.numCodValuesP = StrToInt(tBuf[13]);
                            //std::cout << aBond.atomsCodClasses[0] << std::endl;
                            //std::cout << aBond.atomsCodClasses[1] << std::endl;
                            //std::cout << nline << std::endl;
                            allDictBonds.push_back(aBond);
                         
                            nline+=1;
                            
                            aValueSet tV1;
                            tV1.value    = StrToReal(tBuf[11]);
                            tV1.sigValue = StrToReal(tBuf[12]);
                            tV1.numCodValues = StrToInt(tBuf[13]);
                            allDictBondsIdx1[ha1][ha2][tBuf[2]][tBuf[3]][tBuf[4]][tBuf[5]].push_back(tV1);
                            aValueSet tV2;
                            tV2.value    = StrToReal(tBuf[14]);
                            tV2.sigValue = StrToReal(tBuf[15]);
                            tV2.numCodValues = StrToInt(tBuf[16]);
                            allDictBondsIdx2[ha1][ha2][tBuf[2]][tBuf[3]].push_back(tV2);
                            
                            aValueSet tV3;
                            tV3.value    = StrToReal(tBuf[17]);
                            tV3.sigValue = StrToReal(tBuf[18]);
                            tV3.numCodValues = StrToInt(tBuf[19]);
                            allDictBondsIdx3[ha1][ha2].push_back(tV3);
                           
                        }
                        
                    }
                       
                    codBondFile.close();
                }
                else 
                {
                    std::cout << "Error in setup the programs. " << std::endl;
                    std::cout << iBF->second << " can not be open for reading "
                              << std::endl;
                    exit(1);
                }
                
            }
            catch (std::exception & e)
            {
                std::cout << e.what() << std::endl;
            }
        }
        
        std::cout << "Finish clustering COD org bonds " << std::endl;
        std::time(&tEnd);
        std::cout << "finished at " << std::ctime(&tEnd);
        REAL tDiff;
        tDiff = std::difftime(tEnd,tStart);
        std::cout  << "it takes " << std::setprecision(3) <<tDiff 
                   << " seconds to finish group COD bonds " << std::endl;
        
    }
    
    void  CodClassify::groupCodOrgBonds2()
    {
        setOrgBondHeadHashList2();
        time_t tStart, tEnd;
       
        std::time (&tStart);
        std::cout << "Clustering COD org bonds started at " << std::ctime(&tStart);
        
        int nline = 0;
        for (std::map<ID, ID>::iterator iBF=codOrgBondFiles2.begin();
                    iBF !=codOrgBondFiles2.end(); iBF++)
        {
            try
            {
                // should be something like std::string tNewCodBondFileName(clibMonDir + "/list/bonds.txt");
                // std::string clibMonDir(std::getenv("CLIBD_MON"));
                // std::string codBondFileName = clibMonDir + "/allOrgBonds.table";
                // std::ifstream codBondFile(codBondFileName.c_str());
            
                std::ifstream codBondFile(iBF->second.c_str()); 
                if(codBondFile.is_open())
                {
                    std::string tRecord="";
                    while(!codBondFile.eof())
                    {
                        std::getline(codBondFile, tRecord);
                        
                        tRecord = TrimSpaces(tRecord);
                        std::vector<std::string> tBuf;
                        StrTokenize(tRecord, tBuf);
                        // std::cout << tRecord << std::endl;
                        
                        
                        if ((int)tBuf.size() ==17)
                        {
                            int ha1, ha2;
                        
                            ha1 = StrToInt(tBuf[0]);
                            ha2 = StrToInt(tBuf[1]);
                            
                            
                            allDictBondsIdxD[ha1][ha2][tBuf[2]][tBuf[3]][tBuf[4]][tBuf[5]][tBuf[6]]
                                           [tBuf[7]][tBuf[8]][tBuf[9]][tBuf[10]] = nline;  
                            
                            aValueSet tB;
                            
                            tB.value      = StrToReal(tBuf[11]);
                            tB.sigValue   = StrToReal(tBuf[12]);
                            //aBond.valueP     = StrToReal(tBuf[11]);
                            //aBond.sigValueP  = StrToReal(tBuf[12]);
                            if(tB.sigValue < 0.01)
                            {
                                tB.sigValue = 0.01;
                            }
                            tB.numCodValues  = StrToInt(tBuf[13]);
                            
    
                            allDictBondsD.push_back(tB);
                          
                            
                            nline+=1;
                            //std::cout << tBuf[14] << "  " << tBuf[15] 
                            //          << "  " << tBuf[16] << std::endl;
                            
                           
                            aValueSet tB1;
                            ID a = TrimSpaces(tBuf[14]);
                            tB1.value    = StrToReal(a);
                            tB1.sigValue = StrToReal(tBuf[15]);
                            
                            if (tB1.sigValue <0.01)
                            {
                                tB1.sigValue = 0.01;
                            }
                            tB1.numCodValues = StrToInt(tBuf[16]);
                            
                            allDictBondsIdx1D[ha1][ha2][tBuf[2]][tBuf[3]][tBuf[4]][tBuf[5]]
                                              [tBuf[6]][tBuf[7]][tBuf[8]].push_back(tB1);
                            
                            allDictBondsIdx2D[ha1][ha2][tBuf[2]][tBuf[3]][tBuf[4]][tBuf[5]]
                                              [tBuf[6]].push_back(tB);
                            /*
                            aValueSet tV2;
                            tV2.value    = StrToReal(tBuf[14]);
                            tV2.sigValue = StrToReal(tBuf[15]);
                            tV2.numCodValues = StrToInt(tBuf[16]);
                            allDictBondsIdx2[ha1][ha2][tBuf[2]][tBuf[3]].push_back(tV2);
                            
                            aValueSet tV3;
                            tV3.value    = StrToReal(tBuf[17]);
                            tV3.sigValue = StrToReal(tBuf[18]);
                            tV3.numCodValues = StrToInt(tBuf[19]);
                            allDictBondsIdx3[ha1][ha2].push_back(tV3);
                            */
                           
                        }
                        
                    }
                       
                    codBondFile.close();
                }
                else 
                {
                    std::cout << "Error in setup the programs. " << std::endl;
                    std::cout << iBF->second << " can not be open for reading "
                              << std::endl;
                    exit(1);
                }
                
            }
            catch (std::exception & e)
            {
                std::cout << e.what() << std::endl;
            }
        }
        
        std::cout << "Finish clustering COD org bonds " << std::endl;
        std::time(&tEnd);
        std::cout << "finished at " << std::ctime(&tEnd);
        REAL tDiff;
        tDiff = std::difftime(tEnd,tStart);
        std::cout  << "it takes " << std::setprecision(3) <<tDiff 
                   << " seconds to finish group COD bonds " << std::endl;
        
    }
    
    
    void  CodClassify::groupCodOrgBonds3()
    {
        setOrgBondHeadHashList();
        time_t tStart, tEnd;
       
        std::time (&tStart);
        std::cout << "Clustering COD org bonds started at " << std::ctime(&tStart);
        
        int nline = 0;
        
        for (std::map<ID, ID>::iterator iBF=codOrgBondFiles2.begin();
                    iBF !=codOrgBondFiles2.end(); iBF++)
        {
            try
            {
            
                std::ifstream codBondFile(iBF->second.c_str()); 
                
                if(codBondFile.is_open())
                {
                    std::string tRecord="";
                    while(!codBondFile.eof())
                    {
                        std::getline(codBondFile, tRecord);
                        
                        tRecord = TrimSpaces(tRecord);
                        std::vector<std::string> tBuf;
                        StrTokenize(tRecord, tBuf);
                        
                        if ((int)tBuf.size() ==20)
                        {
                            int ha1, ha2;
                        
                            ha1 = StrToInt(tBuf[0]);
                            ha2 = StrToInt(tBuf[1]);
                            
                            // std::cout << "ha1 " << ha1 << " ha2 " << ha2 << std::endl;
                        
                            /*
                            for (int i=2; i < (int)tBuf.size(); i++)
                            {
                                tBuf[i] = TrimSpaces(tBuf[i]);
                            }
                            */
                            
                            allDictBondsIdx[ha1][ha2][tBuf[2]][tBuf[3]][tBuf[4]][tBuf[5]][tBuf[6]][tBuf[7]]= nline;  
                            
                            BondDict aBond;
                            aBond.seriNum = nline;
                            aBond.atomsHashingCodes.push_back(ha1);
                            aBond.atomsHashingCodes.push_back(ha2);
                            aBond.atomsNB2Rep.push_back(tBuf[2]);
                            aBond.atomsNB2Rep.push_back(tBuf[3]);
                            aBond.atomsNBRep.push_back(tBuf[4]);
                            aBond.atomsNBRep.push_back(tBuf[5]);
                            aBond.atomsCodClasses.push_back(tBuf[6]);
                            aBond.atomsCodClasses.push_back(tBuf[7]);
                            aBond.value      = StrToReal(tBuf[8]);
                            aBond.sigValue   = StrToReal(tBuf[9]);
                            //aBond.valueP     = StrToReal(tBuf[11]);
                            //aBond.sigValueP  = StrToReal(tBuf[12]);
                            if(aBond.sigValue < 0.01)
                            {
                                aBond.sigValue = 0.01;
                            }
                            aBond.numCodValues  = StrToInt(tBuf[10]);
                            //aBond.numCodValuesP = StrToInt(tBuf[13]);
                            //std::cout << aBond.atomsCodClasses[0] << std::endl;
                            //std::cout << aBond.atomsCodClasses[1] << std::endl;
                            //std::cout << nline << std::endl;
                            allDictBonds.push_back(aBond);
                         
                            nline+=1;
                            
                            aValueSet tV1;
                            tV1.value    = StrToReal(tBuf[11]);
                            tV1.sigValue = StrToReal(tBuf[12]);
                            tV1.numCodValues = StrToInt(tBuf[13]);
                            allDictBondsIdx1[ha1][ha2][tBuf[2]][tBuf[3]][tBuf[4]][tBuf[5]].push_back(tV1);
                            aValueSet tV2;
                            tV2.value    = StrToReal(tBuf[14]);
                            tV2.sigValue = StrToReal(tBuf[15]);
                            tV2.numCodValues = StrToInt(tBuf[16]);
                            allDictBondsIdx2[ha1][ha2][tBuf[2]][tBuf[3]].push_back(tV2);
                            
                            aValueSet tV3;
                            tV3.value    = StrToReal(tBuf[17]);
                            tV3.sigValue = StrToReal(tBuf[18]);
                            tV3.numCodValues = StrToInt(tBuf[19]);
                            allDictBondsIdx3[ha1][ha2].push_back(tV3);
                           
                        }
                        
                    }
                       
                    codBondFile.close();
                }
                else 
                {
                    std::cout << "Error in setup the programs. " << std::endl;
                    std::cout << iBF->second << " can not be open for reading "
                              << std::endl;
                    exit(1);
                }
                
            }
            catch (std::exception & e)
            {
                std::cout << e.what() << std::endl;
            }
        }
        
        std::cout << "Finish clustering COD org bonds " << std::endl;
        std::time(&tEnd);
        std::cout << "finished at " << std::ctime(&tEnd);
        REAL tDiff;
        tDiff = std::difftime(tEnd,tStart);
        std::cout  << "it takes " << std::setprecision(3) <<tDiff 
                   << " seconds to finish group COD bonds " << std::endl;
        
    }
    
    void CodClassify::searchCodBonds()
    {
        int nFind =0;
              
        for (std::vector<BondDict>::iterator iGB = allBonds.begin();
                iGB != allBonds.end(); iGB++)
        {
            std::cout << "Bond between "<< iGB->atoms[0] << " and "
                    << iGB->atoms[1] << std::endl; 
            
            bool mB = false;
            for (std::map<ID, int>::iterator iKey=iGB->fullAtoms.begin();
                    iKey !=iGB->fullAtoms.end(); iKey++)
            {
                if (allAtoms[iKey->second].isMetal)
                {
                    mB = true;
                }
            }
            
            
            
            if (mB)
            {
                std::cout << " search a metal-related bond " << std::endl;
                exit(1);
                
                searchCodMetBonds(iGB);
            }
            else
            {
                std::cout << " search a organic-only bond " << std::endl;
                
                searchCodOrgBonds2(iGB);
            }
            
            std::cout << "The final target bond value is " 
                      << iGB->value << std::endl;
            if (iGB->value > 0.3)
            {
                nFind++;
            }
            
            
        }   
        std::cout << "Number of Bonds to be found " << (int)allBonds.size() << std::endl
                  << "Number of Bonds found " << nFind << std::endl;
        
    }
    
    void CodClassify::searchCodOrgBonds(std::vector<BondDict>::iterator iB)
    {   
        
            std::vector<int> tPair;
            for (std::map<ID, int>::iterator iA=iB->fullAtoms.begin();
                    iA !=iB->fullAtoms.end(); iA++)
            {
                tPair.push_back(iA->second);
            }
            int ha1, ha2;
            ID a1NB2, a2NB2, a1NB, a2NB, a1C, a2C;
            
            if((int) allAtoms[tPair[0]].hashingValue < 
                   (int) allAtoms[tPair[1]].hashingValue )
            {
                ha1  =allAtoms[tPair[0]].hashingValue;
                ha2  =allAtoms[tPair[1]].hashingValue;
                a1NB2= allAtoms[tPair[0]].codNB2Symb;
                a2NB2= allAtoms[tPair[1]].codNB2Symb;
                a1NB = allAtoms[tPair[0]].codNBSymb;
                a2NB = allAtoms[tPair[1]].codNBSymb;
                a1C  = allAtoms[tPair[0]].codClass;
                a2C  = allAtoms[tPair[1]].codClass;
            }
            else if ((int) allAtoms[tPair[0]].hashingValue == 
                   (int) allAtoms[tPair[1]].hashingValue)
            {
                if((int)allAtoms[tPair[0]].codClass.size() <= 
                   (int)allAtoms[tPair[1]].codClass.size())
                {
                    ha1  =allAtoms[tPair[0]].hashingValue;
                    ha2  =allAtoms[tPair[1]].hashingValue;
                    a1NB2= allAtoms[tPair[0]].codNB2Symb;
                    a2NB2= allAtoms[tPair[1]].codNB2Symb;
                    a1NB = allAtoms[tPair[0]].codNBSymb;
                    a2NB = allAtoms[tPair[1]].codNBSymb;
                    a1C  = allAtoms[tPair[0]].codClass;
                    a2C  = allAtoms[tPair[1]].codClass;
                }
                else
                {
                    ha1  =allAtoms[tPair[1]].hashingValue;
                    ha2  =allAtoms[tPair[0]].hashingValue;
                    a1NB2= allAtoms[tPair[1]].codNB2Symb;
                    a2NB2= allAtoms[tPair[0]].codNB2Symb;
                    a1NB = allAtoms[tPair[1]].codNBSymb;
                    a2NB = allAtoms[tPair[0]].codNBSymb;
                    a1C  = allAtoms[tPair[1]].codClass;
                    a2C  = allAtoms[tPair[0]].codClass;
                }
            }
            else
            {
                ha1  =allAtoms[tPair[1]].hashingValue;
                ha2  =allAtoms[tPair[0]].hashingValue;
                a1NB2= allAtoms[tPair[1]].codNB2Symb;
                a2NB2= allAtoms[tPair[0]].codNB2Symb;
                a1NB = allAtoms[tPair[1]].codNBSymb;
                a2NB = allAtoms[tPair[0]].codNBSymb;
                a1C  = allAtoms[tPair[1]].codClass;
                a2C  = allAtoms[tPair[0]].codClass;
            }
            
            int dLev = 0;
            
            std::cout << "for target bond of atoms " <<iB->atoms[0] << " and "
                        << iB->atoms[1] <<std::endl;
                       
            for (std::map<ID, int>::iterator iAt=iB->fullAtoms.begin();
                        iAt != iB->fullAtoms.end(); iAt++)
            {
                std::cout << "chemType : " << allAtoms[iAt->second].chemType << " COD classes : " 
                          << allAtoms[iAt->second].codClass << std::endl;
            }
                
            std::cout << "ha1 " << ha1 << " ha2 " << ha2 << std::endl
                      << " a1NB2 " << a1NB2 << " a2NB2 " << a2NB2  << std::endl
                      << " a1NB "  << a1NB  << " a2NB " << a2NB << std::endl
                      << " a1C "   << a1C   << " a2C "  << a2C << std::endl;
            
            
            if ((int)allDictBondsIdx[ha1][ha2].size() !=0)
            {
                std::cout << "Found all hashing codes" << std::endl;
              
                std::map<ID, std::map<ID,  std::map<ID,
                std::map<ID, std::map<ID, std::map<ID,  
                int> > > > > >::iterator iFind1 = allDictBondsIdx[ha1][ha2].find(a1NB2);
                if (iFind1 !=allDictBondsIdx[ha1][ha2].end())
                {
                    // std::cout << "find " << a1NB2 << std::endl;
                    
                    std::map<ID, std::map<ID,  std::map<ID, std::map<ID, std::map<ID, 
                    int> > > > >::iterator iFind2 = allDictBondsIdx[ha1][ha2][a1NB2].find(a2NB2);
                    
                    // std::cout << "size " << (int)allDictBondsIdx[ha1][ha2][a1NB2].size() << std::endl;
                    
                    if (iFind2 != allDictBondsIdx[ha1][ha2][a1NB2].end())
                    {
                       //std::cout << " find " << a2NB2 << std::endl; 
                       // std::cout << "size " << (int)allDictBondsIdx[ha1][ha2][a1NB2][a2NB2].size() << std::endl;
                        
                        std::map<ID, std::map<ID,  std::map<ID, std::map<ID,
                        int> > > >::iterator iFind3 = allDictBondsIdx[ha1][ha2][a1NB2][a2NB2].find(a1NB);
                        if (iFind3 != allDictBondsIdx[ha1][ha2][a1NB2][a2NB2].end())
                        {
                            //std::cout << " find " << a1NB << std::endl; 
                            // std::cout << "size " << (int)allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB].size() << std::endl;
                            
                            std::map<ID, std::map<ID,  std::map<ID, int> > >::iterator iFind4
                                = allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB].find(a2NB);
                            if(iFind4 !=allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB].end())
                            {
                               /*
                               std::cout << " find " << a2NB << std::endl; 
                               std::cout << "size " << (int)allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB].size() 
                                        << std::endl;
                                */
                                std::map<ID, std::map<ID,int> >::iterator iFind5
                                        = allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB].find(a1C);
                                if (iFind5 != allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB].end())
                                {
                                    /*
                                    std::cout << " find " << a1C << std::endl; 
                                    std::cout << "size " << (int)allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][a1C].size() 
                                              << std::endl;
                                     */
                                    std::map<ID, int>::iterator iFind6 
                                    = allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][a1C].find(a2C);
                                    if (iFind6 !=allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][a1C].end())
                                    {
                                        
                                        
                                        
                                        int tIdx = allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][a1C][a2C];
                                        /*
                                        std::cout << " find " << a2C << std::endl; 
                                        std::cout << "size "  << allDictBonds[tIdx].numCodValues << std::endl;
                                         */
                                        // find a COD bond with the exact class 
                                        std::cout << "Found exact matches of atom cod-classes " << std::endl;
                                        if (allDictBonds[tIdx].numCodValues >= 5)
                                        {
                                            iB->hasCodValue = true;
                                            iB->value = allDictBonds[tIdx].value;
                                            iB->sigValue = allDictBonds[tIdx].sigValue;
                                        }
                                        else
                                        {
                                         
                                            std::vector<BondDict> tBs6;
                                            for (std::map<ID, std::map<ID, int> >::iterator iB5
                                            =allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB].begin();
                                            iB5 !=allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB].end();
                                            iB5++)
                                            {
                                                for (std::map<ID, int>::iterator iB6=iB5->second.begin();
                                                     iB6 !=iB5->second.end(); iB6++)
                                                {
                                                    tBs6.push_back(allDictBonds[iB6->second]);
                                                }
                                            }
                                            setupTargetBondsUsingMean(tBs6, iB);
                                        }
                                    }
                                    else  // one atom without exact matching of the cod class, iFind6.
                                    {
                                        std::vector<BondDict> tBs6;
                                        for (std::map<ID, int>::iterator iB6
                                                =allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][a1C].begin();
                                                iB6 !=allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][a1C].end();
                                                iB6++)
                                        {
                                            tBs6.push_back(allDictBonds[iB6->second]);
                                        }
                                        // dissimilarity is small, using mean values on this cluster of atoms
                                        setupTargetBondsUsingMean(tBs6, iB);
                                    }
                                }
                                else  // both atoms without exact matching of cod classes, iFind5.
                                {
                                    std::vector<BondDict> tBs6;
                                    for (std::map<ID, std::map<ID, int> >::iterator iB5
                                          =allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB].begin();
                                         iB5 !=allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB].end();
                                         iB5++)
                                    {
                                        for (std::map<ID, int>::iterator iB6=iB5->second.begin();
                                                iB6 !=iB5->second.end(); iB6++)
                                        {
                                            tBs6.push_back(allDictBonds[iB6->second]);
                                        }
                                    }
                                    // dissimilarity is still small, using mean values on this cluster of atoms
                                    setupTargetBondsUsingMean(tBs6, iB);
                                    
                                }
                            }
                            else // one 2NB space not matching (both 2NB conf not match) iFind4 
                            {
                                
                                std::vector<BondDict> tBs6;
                                
                                //std::cout << (int)allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB].size()
                                //       << std::endl;
                                
                                for (std::map<ID, std::map<ID, std::map<ID, int> > >::iterator iB4
                                          =allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB].begin();
                                         iB4 !=allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB].end();
                                         iB4++)
                                {
                                    
                                    for(std::map<ID, std::map<ID, int> >::iterator iB5=iB4->second.begin();
                                            iB5!=iB4->second.end(); iB5++)
                                    {
                                        for (std::map<ID, int>::iterator iB6=iB5->second.begin();
                                                iB6 !=iB5->second.end(); iB6++)
                                        {
                                            tBs6.push_back(allDictBonds[iB6->second]);
                                            
                                        }
                                    }
                                }
                               
                                // one atom pair does not match in 2NB conf
                                dLev =1; 
                                
                                setupTargetBondsUsingSymblDist(tBs6, iB, dLev);
                            }
                        }
                        else // both atom pair do not match 2NB conf not matching, iFind3 
                        {
                            
                            std::vector<BondDict> tBs6;
                            for (std::map<ID, std::map<ID, std::map<ID, std::map<ID, int> > > >::iterator iB3
                                 =allDictBondsIdx[ha1][ha2][a1NB2][a2NB2].begin();
                                 iB3 !=allDictBondsIdx[ha1][ha2][a1NB2][a2NB2].end();
                                 iB3++)
                            {
                                for(std::map<ID, std::map<ID, std::map<ID, int> > >::iterator iB4=iB3->second.begin();
                                    iB4!=iB3->second.end(); iB4++)
                                {
                                    for(std::map<ID, std::map<ID, int> >::iterator iB5=iB4->second.begin();
                                        iB5!=iB4->second.end(); iB5++)
                                    {
                                        for (std::map<ID, int>::iterator iB6=iB5->second.begin();
                                                iB6 !=iB5->second.end(); iB6++)
                                        {
                                            tBs6.push_back(allDictBonds[iB6->second]);
                                        }
                                    }
                                }
                            }
                            dLev = 1;
                            setupTargetBondsUsingSymblDist(tBs6, iB, dLev);
                        }
                    }
                    else // One atom pair does not match 2NB space (both not match 2NB confs), iFind2
                    {
                        std::vector<BondDict> tBs6;
                        for (std::map<ID, std::map<ID, std::map<ID, std::map<ID, std::map<ID, int> > > > >::iterator iB2
                                 =allDictBondsIdx[ha1][ha2][a1NB2].begin();
                                 iB2 !=allDictBondsIdx[ha1][ha2][a1NB2].end();
                                 iB2++)
                        {
                            for(std::map<ID, std::map<ID, std::map<ID, std::map<ID, int> > > >::iterator iB3=iB2->second.begin();
                                    iB3!=iB2->second.end(); iB3++)
                            {
                                for(std::map<ID, std::map<ID, std::map<ID, int> > >::iterator iB4=iB3->second.begin();
                                    iB4!=iB3->second.end(); iB4++)
                                {
                                    for(std::map<ID, std::map<ID, int> >::iterator iB5=iB4->second.begin();
                                        iB5!=iB4->second.end(); iB5++)
                                    {
                                        for (std::map<ID, int>::iterator iB6=iB5->second.begin();
                                                iB6 !=iB5->second.end(); iB6++)
                                        {
                                            tBs6.push_back(allDictBonds[iB6->second]);
                                        }
                                    }
                                }
                            }
                        }
                        dLev = 2;
                        setupTargetBondsUsingSymblDist(tBs6, iB, dLev);
                    }
                }
                else // Both atom pairs do not match 2NB spaces, corresponding to iFind1
                {
                    std::vector<BondDict> tBs6;
                    for (std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                            std::map<ID, std::map<ID, int> > > > > >::iterator iB1
                            =allDictBondsIdx[ha1][ha2].begin();
                            iB1 !=allDictBondsIdx[ha1][ha2].end();
                            iB1++)
                    {
                        for(std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                                std::map<ID, int> > > > >::iterator iB2=iB1->second.begin();
                                    iB2!=iB1->second.end(); iB2++)
                        {
                            for(std::map<ID, std::map<ID, std::map<ID, std::map<ID, int> > > >::iterator iB3=iB2->second.begin();
                                    iB3!=iB2->second.end(); iB3++)
                            {
                                for(std::map<ID, std::map<ID, std::map<ID, int> > >::iterator iB4=iB3->second.begin();
                                    iB4!=iB3->second.end(); iB4++)
                                {
                                    for(std::map<ID, std::map<ID, int> >::iterator iB5=iB4->second.begin();
                                        iB5!=iB4->second.end(); iB5++)
                                    {
                                        for (std::map<ID, int>::iterator iB6=iB5->second.begin();
                                                iB6 !=iB5->second.end(); iB6++)
                                        {
                                            tBs6.push_back(allDictBonds[iB6->second]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    dLev = 2;
                    setupTargetBondsUsingSymblDist(tBs6, iB, dLev);
                }
            }
            else
            {
                std::cout << "Could not find the basic bond classes in Cod " << std::endl
                        << "for dictionary bond of atoms " <<iB->atoms[0] << " and "
                        << iB->atoms[1] <<std::endl;
                
                for (std::map<ID, int>::iterator iAt=iB->fullAtoms.begin();
                        iAt != iB->fullAtoms.end(); iAt++)
                {
                    std::cout << "chemType : " << allAtoms[iAt->second].chemType << " atom COD classes : " 
                              << allAtoms[iAt->second].codClass << std::endl;
                }
            }
            
            
            //std::cout << "The bond value is now " << iB->value << std::endl;
            //if (iB->hasCodValue)
            //{
            //    std::cout << "It have exact matches to COD atoms " << std::endl;
            //}
            //else
            //{
            //    std::cout << "The bond value is generated by " << std::endl
            //              << "either averaging or shortest distance to COD"
            //              << std::endl;
            //}
    }
    
    void CodClassify::searchCodOrgBonds2(std::vector<BondDict>::iterator iB)
    {   
            std::vector<int> tPair;
            for (std::map<ID, int>::iterator iA=iB->fullAtoms.begin();
                    iA !=iB->fullAtoms.end(); iA++)
            {
                tPair.push_back(iA->second);
            }
            
            int ha1, ha2;
            ID a1NB2, a2NB2, a1NB, a2NB, a1M, a2M, a1C, a2C;
            // int as0=-1, as1=-1;
            
            if((int) allAtoms[tPair[0]].hashingValue < 
                   (int) allAtoms[tPair[1]].hashingValue )
            {
                ha1  =allAtoms[tPair[0]].hashingValue;
                ha2  =allAtoms[tPair[1]].hashingValue;
                a1NB2= allAtoms[tPair[0]].codNB2Symb;
                a2NB2= allAtoms[tPair[1]].codNB2Symb;
                a1NB = allAtoms[tPair[0]].codNBSymb;
                a2NB = allAtoms[tPair[1]].codNBSymb;
                a1M  = allAtoms[tPair[0]].codAtmMain;
                a2M  = allAtoms[tPair[1]].codAtmMain;
                a1C  = allAtoms[tPair[0]].codClass;
                a2C  = allAtoms[tPair[1]].codClass;
                //as0 =0;
                //as1 =1;
            }
            else if ((int) allAtoms[tPair[0]].hashingValue == 
                   (int) allAtoms[tPair[1]].hashingValue)
            {
                std::vector<sortMap4> vM;
                
                sortMap4   m1, m2;
                
                m1.id    = allAtoms[tPair[0]].id;
                m2.id    = allAtoms[tPair[0]].id;
                m1.ha    =(int)allAtoms[tPair[0]].hashingValue;
                m2.ha    =(int)allAtoms[tPair[1]].hashingValue;
                m1.lev2  = allAtoms[tPair[0]].codNB2Symb;
                m2.lev2  = allAtoms[tPair[1]].codNB2Symb;
                m1.lev3  = allAtoms[tPair[0]].codNBSymb;
                m2.lev3  = allAtoms[tPair[1]].codNBSymb;
                m1.key   = allAtoms[tPair[0]].codAtmMain;
                m2.key   = allAtoms[tPair[1]].codAtmMain;
                m1.lev4  = allAtoms[tPair[0]].codClass;
                m2.lev4  = allAtoms[tPair[1]].codClass;
                
                vM.push_back(m1);
                vM.push_back(m2);
                
                std::sort(vM.begin(), vM.end(), sortMapkey4);
                
                ha1    = vM[0].ha;
                ha2    = vM[1].ha;
                a1NB2  = vM[0].lev2;
                a2NB2  = vM[1].lev2;
                a1NB   = vM[0].lev3;
                a2NB   = vM[1].lev3;
                a1M    = vM[0].key;
                a2M    = vM[1].key;
                a1C    = vM[0].lev4;
                a2C    = vM[1].lev4;
                
                //as0 =0;
                //as1 =1;
                
                
                /*
                if((int)allAtoms[tPair[0]].codClass.size() <= 
                   (int)allAtoms[tPair[1]].codClass.size())
                {
                    ha1  =allAtoms[tPair[0]].hashingValue;
                    ha2  =allAtoms[tPair[1]].hashingValue;
                    a1NB2= allAtoms[tPair[0]].codNB2Symb;
                    a2NB2= allAtoms[tPair[1]].codNB2Symb;
                    a1NB = allAtoms[tPair[0]].codNBSymb;
                    a2NB = allAtoms[tPair[1]].codNBSymb;
                    a1C  = allAtoms[tPair[0]].codClass;
                    a2C  = allAtoms[tPair[1]].codClass;
                    as0=0;
                    as1=1;
                }
                else
                {
                    ha1  =allAtoms[tPair[1]].hashingValue;
                    ha2  =allAtoms[tPair[0]].hashingValue;
                    a1NB2= allAtoms[tPair[1]].codNB2Symb;
                    a2NB2= allAtoms[tPair[0]].codNB2Symb;
                    a1NB = allAtoms[tPair[1]].codNBSymb;
                    a2NB = allAtoms[tPair[0]].codNBSymb;
                    a1C  = allAtoms[tPair[1]].codClass;
                    a2C  = allAtoms[tPair[0]].codClass;
                    as0 =1;
                    as1 =0;
                }
                 */
            }
            else
            {
                ha1  =allAtoms[tPair[1]].hashingValue;
                ha2  =allAtoms[tPair[0]].hashingValue;
                a1NB2= allAtoms[tPair[1]].codNB2Symb;
                a2NB2= allAtoms[tPair[0]].codNB2Symb;
                a1NB = allAtoms[tPair[1]].codNBSymb;
                a2NB = allAtoms[tPair[0]].codNBSymb;
                a1M  = allAtoms[tPair[1]].codAtmMain;
                a2M  = allAtoms[tPair[0]].codAtmMain;
                a1C  = allAtoms[tPair[1]].codClass;
                a2C  = allAtoms[tPair[0]].codClass;
                //as0 =1;
                //as1 =0;
            }
            
            ID tInR, tInR1("Y"), tInR2("N");
            
            if (iB->isInSameRing)
            {
                tInR = "Y";
            }
            else
            {
                tInR = "N";
            }
            
            // int dLev = 0;
            
            std::cout << "for target bond of atoms " <<iB->atoms[0] << " and "
                        << iB->atoms[1] <<std::endl;
                       
            for (std::map<ID, int>::iterator iAt=iB->fullAtoms.begin();
                        iAt != iB->fullAtoms.end(); iAt++)
            {
                std::cout << "chemType : " << allAtoms[iAt->second].chemType << " COD classes : " 
                          << allAtoms[iAt->second].codClass << std::endl;
            }
                
            std::cout << "ha1 " << ha1 << " ha2 " << ha2 << std::endl
                      << " a1NB2 " << a1NB2 << " a2NB2 " << a2NB2  << std::endl
                      << " a1NB "  << a1NB  << " a2NB " << a2NB << std::endl
                      << " a1M " << a1M     << " a2M "  << a2M <<  std::endl
                      << " a1C "   << a1C   << " a2C "  << a2C << std::endl;
            
            std::cout << "atom 1 " << allAtoms[tPair[0]].ccp4Type
                      << "  atom 2 " << allAtoms[tPair[1]].ccp4Type
                      << std::endl;
            
            if (allDictBondsIdxD.find(ha1) != allDictBondsIdxD.end()
                && allDictBondsIdxD[ha1].find(ha2) != allDictBondsIdxD[ha1].end())
            {
                std::cout << "Found all hashing codes " << std::endl;
                
                // std::map<ID, std::map<ID,
                //     std::vector<aValueSet> > >::iterator iFind1 = allDictBondsIdx2[ha1][ha2].find(a1NB2);
                if ( allDictBondsIdxD[ha1][ha2].find(a1NB2) !=allDictBondsIdxD[ha1][ha2].end())
                {
                    std::cout << "find " << a1NB2 << std::endl;
                    // iFind1 
                    // std::map<ID,  std::vector<aValueSet> >::iterator iFind2 = allDictBondsIdx2[ha1][ha2][a1NB2].find(a2NB2);
                    
                    //std::cout << "size " << (int)allDictBondsIdx[ha1][ha2][a1NB2].size() << std::endl;
                    
                    if (allDictBondsIdxD[ha1][ha2][a1NB2].find(a2NB2) != allDictBondsIdxD[ha1][ha2][a1NB2].end())
                    {
                        std::cout << " find " << a2NB2 << std::endl; 
                        // std::cout << "size " << (int)allDictBondsIdx[ha1][ha2][a1NB2][a2NB2].size() << std::endl;
                        // iFind2  
                        //std::map<ID, std::map<ID,  std::map<ID, std::map<ID,
                        //int> > > >::iterator iFind3 = allDictBondsIdx[ha1][ha2][a1NB2][a2NB2].find(a1NB);
                        if ( allDictBondsIdxD[ha1][ha2][a1NB2][a2NB2].find(a1NB)
                                != allDictBondsIdxD[ha1][ha2][a1NB2][a2NB2].end())
                        {
                            // iFind3
                            std::cout << " find " << a1NB << std::endl; 
                            // std::cout << "size " << (int)allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB].size() << std::endl;
                            
                            //std::map<ID, std::map<ID,  std::map<ID, int> > >::iterator iFind4
                            //    = allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB].find(a2NB);
                            if( allDictBondsIdxD[ha1][ha2][a1NB2][a2NB2][a1NB].find(a2NB)
                                 !=allDictBondsIdxD[ha1][ha2][a1NB2][a2NB2][a1NB].end())
                            {
                                // iFind4
                                std::cout << " find " << a2NB << std::endl;
                               
                                //std::cout << "size " << (int)allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB].size() 
                                //         << std::endl;
                                
                                // std::map<ID, std::map<ID,int> >::iterator iFind5
                                //         = allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB].find(a1C)find(a1C);
                                
                                if (allDictBondsIdxD[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB].find(tInR1)
                                         != allDictBondsIdxD[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB].end())
                                {
                                    tInR = tInR1;
                                }
                                else if (allDictBondsIdxD[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB].find(tInR2)
                                         != allDictBondsIdxD[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB].end())
                                {
                                    tInR = tInR2;
                                }
                                else
                                {
                                    std::cout << "Acedrg database error for the bond between atoms "
                                              << a1C << " and " << a2C << std::endl;
                                    exit(1);
                                }
                               
                                if ( allDictBondsIdxD[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][tInR].find(a1M)
                                         != allDictBondsIdxD[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][tInR].end())
                                {
                                    // iFind 5 a1M
                                    std::cout << " find " << a1M << std::endl; 
                                    //std::cout << "size " << (int)allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][a1C].size() 
                                    //          << std::endl;
                                    
                                    // std::map<ID, int>::iterator iFind6 a2M 
                                    // = allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][a1C].find(a2C);
                                    if( allDictBondsIdxD[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][tInR][a1M].find(a2M)
                                         !=allDictBondsIdxD[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][tInR][a1M].end())
                                    {
                                        // iFind6
                                        std::cout << " find " << a2M << std::endl;
                                        if (allDictBondsIdxD[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][tInR][a1M][a2M].find(a1C)
                                            !=allDictBondsIdxD[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][tInR][a1M][a2M].end())
                                        {
                                            //iFind7 a1C
                                            std::cout << " find " << a1C << std::endl;
                                            if (allDictBondsIdxD[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][tInR][a1M][a2M][a1C].find(a2C)
                                            !=allDictBondsIdxD[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][tInR][a1M][a2M][a1C].end())
                                            {
                                                // iFind 8 a2C
                                                std::cout << " find " << a2C << std::endl;
                                                int tIdx = allDictBondsIdxD[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][tInR][a1M][a2M][a1C][a2C];
                                                /*
                                                   std::cout << " find " << a2C << std::endl; 
                                                   std::cout << "size "  << allDictBonds[tIdx].numCodValues << std::endl;
                                                */
                                                // find a COD bond with the exact class 
                                                std::cout << "Found exact matches of atom cod-classes " << std::endl;
                                                if (allDictBondsD[tIdx].numCodValues > 5)
                                                {
                                                    iB->hasCodValue = true;
                                                    iB->value    = allDictBondsD[tIdx].value;
                                                    iB->sigValue = allDictBondsD[tIdx].sigValue;  
                                                }
                                                else
                                                {
                                                    iB->hasCodValue = true;
                                                    iB->value    = allDictBondsIdx1D[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][tInR][a1M][a2M][0].value;
                                                    iB->sigValue = allDictBondsIdx1D[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][tInR][a1M][a2M][0].sigValue;
                                                }
                                            }
                                            else
                                            {
                                                // iFind 8 failed a2C
                                                iB->hasCodValue = true;
                                                iB->value    = allDictBondsIdx1D[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][tInR][a1M][a2M][0].value;
                                                iB->valueST  =iB->value;
                                                iB->sigValue = allDictBondsIdx1D[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][tInR][a1M][a2M][0].sigValue;
                                                iB->sigValueST =iB->sigValue;
                                                std::cout << "iFind 8" << std::endl;
                                            }
                                            
                                        }
                                        else
                                        {
                                            // iFind 7 failed a1C
                                            iB->hasCodValue = true;
                                            iB->value    = allDictBondsIdx1D[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][tInR][a1M][a2M][0].value;
                                            iB->valueST  =iB->value;
                                            iB->sigValue = allDictBondsIdx1D[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][tInR][a1M][a2M][0].sigValue;
                                            iB->valueST  =iB->value;
                                            std::cout << "iFind 7" << std::endl;
                                        }        
                                    }
                                    else   
                                    {
                                        // iFind6 failed a2M
                                        aValueSet   tVaS;
                                        setValueSet(tVaS, allDictBondsIdx2D[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][tInR]);
                                        iB->value        =tVaS.value;
                                        iB->valueST      =iB->value;
                                        iB->sigValue     =tVaS.sigValue;
                                        iB->sigValueST   =iB->sigValue;
                                        iB->numCodValues =tVaS.numCodValues;
                                        std::cout << "iFind 6" << std::endl;
                                    }
                                }
                                else  // iFind5 failed a1M .
                                {
                                    aValueSet   tVaS;
                                    setValueSet(tVaS, allDictBondsIdx2D[ha1][ha2][a1NB2][a2NB2][a1NB][a2NB][tInR]);
                                    iB->value        =tVaS.value;
                                    iB->valueST      =iB->value;
                                    iB->sigValue     =tVaS.sigValue;
                                    iB->sigValueST   =iB->sigValue;
                                    iB->numCodValues =tVaS.numCodValues;
                                    std::cout << "iFind 5: " << a1M << std::endl;
                                }
                            }
                            else // iFind4 failed  a2NB
                            {
                                
                                std::vector<aValueSet> tBs4;
                                
                                //std::cout << (int)allDictBondsIdx[ha1][ha2][a1NB2][a2NB2][a1NB].size()
                                //       << std::endl;
                                
                                for (std::map<ID, std::map<ID, std::vector<aValueSet> > >::iterator iB4
                                          =allDictBondsIdx2D[ha1][ha2][a1NB2][a2NB2][a1NB].begin();
                                         iB4 !=allDictBondsIdx2D[ha1][ha2][a1NB2][a2NB2][a1NB].end();
                                         iB4++)
                                {
                                    for (std::map<ID, std::vector<aValueSet> >::iterator iB5=iB4->second.begin();
                                            iB5 !=iB4->second.end(); iB5++)
                                    {
                                        for(std::vector<aValueSet>::iterator iB6=iB5->second.begin();
                                                iB6 !=iB5->second.end(); iB6++)
                                        {
                                            tBs4.push_back(*iB6);
                                        }
                                    }
                                }
                                
                                aValueSet   tVaS;
                                setValueSet(tVaS, tBs4);
                                iB->value        =tVaS.value;
                                iB->valueST      =iB->value;
                                iB->sigValue     =tVaS.sigValue;
                                iB->sigValueST   =iB->sigValue;
                                iB->numCodValues =tVaS.numCodValues;                
                                std::cout << "iFind 4" << std::endl;
                            }
                        }
                        else // iFind3 failed a1NB2 
                        {
                            std::vector<aValueSet> tBs3;
                                
                            for (std::map<ID, std::map<ID, std::map<ID, std::vector<aValueSet> > > >::iterator iB3
                                     =allDictBondsIdx2D[ha1][ha2][a1NB2][a2NB2].begin();
                                 iB3 !=allDictBondsIdx2D[ha1][ha2][a1NB2][a2NB2].end();
                                         iB3++)
                            {
                                for (std::map<ID, std::map<ID, std::vector<aValueSet> > >::iterator iB4
                                        =iB3->second.begin(); iB4 !=iB3->second.end(); iB4++)
                                {
                                    for (std::map<ID, std::vector<aValueSet> >::iterator iB5=iB4->second.begin();
                                        iB5 !=iB4->second.end(); iB5++)
                                    {
                                        for(std::vector<aValueSet>::iterator iB6=iB5->second.begin();
                                                iB6 !=iB5->second.end(); iB6++)
                                        {
                                            tBs3.push_back(*iB6);
                                        }
                                    }
                                }
                            }
                            
                            aValueSet   tVaS;
                            setValueSet(tVaS, tBs3);
                            
                            std::cout << "iFind  3 A1 " << std::endl
                                      << "tVaS num " << tVaS.numCodValues << std::endl
                                      << "tVaS sig " << tVaS.sigValue << std::endl;
                            
                            if (tVaS.numCodValues > 10 && tVaS.sigValue <0.05)
                            {
                                iB->value        =tVaS.value;
                                iB->valueST      =iB->value;
                                iB->sigValue     =tVaS.sigValue;
                                iB->sigValueST   =iB->sigValue;
                                iB->numCodValues =tVaS.numCodValues;
                                
                                std::cout << "iFind 3 A" << std::endl;
                                
                            }
                            else
                            {
                                //Test 
                                for (std::map<ID, std::map<ID, std::map<ID, 
                                 std::map<ID, std::vector<aValueSet> > > > >::iterator iB3
                                     =allDictBondsIdx2D[ha1][ha2][a1NB2].begin();
                                     iB3 !=allDictBondsIdx2D[ha1][ha2][a1NB2].end();
                                         iB3++)
                                {
                                    if (iB3->first !=a2NB2)
                                    {
                                        for (std::map<ID, std::map<ID, std::map<ID,  
                                             std::vector<aValueSet> > > >::iterator iB4
                                             =iB3->second.begin(); iB4 !=iB3->second.end(); iB4++)
                                        {
                                            for (std::map<ID, std::map<ID, std::vector<aValueSet> > >::iterator iB5=iB4->second.begin();
                                                 iB5 !=iB4->second.end(); iB5++)
                                            {
                                                for(std::map<ID, std::vector<aValueSet> >::iterator iB6=iB5->second.begin();
                                                    iB6 !=iB5->second.end(); iB6++)
                                                {
                                                    for(std::vector<aValueSet>::iterator iB7=iB6->second.begin();
                                                        iB7 !=iB6->second.end(); iB7++)
                                                    {
                                                        tBs3.push_back(*iB7);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                
                                for (std::map<ID, std::map<ID, std::map<ID, 
                                     std::map<ID, std::map<ID, std::vector<aValueSet> > > > > >::iterator iB1
                                     =allDictBondsIdx2D[ha1][ha2].begin();
                                     iB1 !=allDictBondsIdx2D[ha1][ha2].end(); iB1++)
                                {
                                    if (iB1->first !=a1NB2)
                                    {
                                        for (std::map<ID, std::map<ID, std::map<ID, 
                                             std::map<ID, std::vector<aValueSet> > > > >::iterator iB2
                                             =iB1->second.begin(); 
                                             iB2 !=iB1->second.end(); iB2++)
                                        {
                                            if (iB2->first == a2NB2)
                                            {
                                                for (std::map<ID, std::map<ID, std::map<ID,  
                                                     std::vector<aValueSet> > > >::iterator iB3
                                                     =iB2->second.begin(); iB3 !=iB2->second.end(); iB3++)
                                                {
                                                    for (std::map<ID, std::map<ID, std::vector<aValueSet> > >::iterator iB4=iB3->second.begin();
                                                         iB4 !=iB3->second.end(); iB4++)
                                                    {
                                                        for(std::map<ID, std::vector<aValueSet> >::iterator iB5=iB4->second.begin();
                                                            iB5 !=iB4->second.end(); iB5++)
                                                        {
                                                            for(std::vector<aValueSet>::iterator iB6=iB5->second.begin();
                                                                iB6 !=iB5->second.end(); iB6++)
                                                            {
                                                                tBs3.push_back(*iB6);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                
                                setValueSet(tVaS, tBs3);
                                std::cout << "iFind  3 A2 " << std::endl
                                        << "tVaS num " << tVaS.numCodValues << std::endl
                                        << "tVaS sig " << tVaS.sigValue << std::endl;
                                
                                        
                                if (tVaS.numCodValues > 10 && tVaS.sigValue <0.05)
                                {
                                    iB->value        =tVaS.value;
                                    iB->valueST      =iB->value;
                                    iB->sigValue     =tVaS.sigValue;
                                    iB->sigValueST   =iB->sigValue;
                                    iB->numCodValues =tVaS.numCodValues;
                                
                                    std::cout << "iFind 3 A2" << std::endl;
                                
                                }
                                else if (ccp4BondsA.find(allAtoms[tPair[0]].ccp4Type) !=ccp4BondsA.end() 
                                         && ccp4BondsA[allAtoms[tPair[0]].ccp4Type].find(allAtoms[tPair[1]].ccp4Type)
                                         !=ccp4BondsA[allAtoms[tPair[0]].ccp4Type].end())
                                {
                                    getCCP4Bonds(iB, allAtoms[tPair[0]].ccp4Type, allAtoms[tPair[1]].ccp4Type);
                                    std::cout << "iFind 3 B" << std::endl;
                                }
                                else if (ccp4BondsA.find(allAtoms[tPair[1]].ccp4Type) !=ccp4BondsA.end() 
                                         && ccp4BondsA[allAtoms[tPair[1]].ccp4Type].find(allAtoms[tPair[0]].ccp4Type)
                                         !=ccp4BondsA[allAtoms[tPair[1]].ccp4Type].end())
                                {
                                    getCCP4Bonds(iB, allAtoms[tPair[1]].ccp4Type, allAtoms[tPair[0]].ccp4Type);
                                    std::cout << "iFind 3 C" << std::endl;
                                }
                                else
                                {
                                    iB->value        =tVaS.value;
                                    iB->valueST      =iB->value;
                                    iB->sigValue     =tVaS.sigValue;
                                    iB->sigValueST   =iB->sigValue;
                                    iB->numCodValues =tVaS.numCodValues;   
                                    std::cout << "iFind 3" << std::endl;    
                                }
                            }
                        }
                    }
                    else // iFind2 a2NB2 
                    {
                        std::vector<aValueSet> tBs2;
                                
                        for (std::map<ID, std::map<ID, std::map<ID, 
                             std::map<ID, std::vector<aValueSet> > > > >::iterator iB3
                                     =allDictBondsIdx2D[ha1][ha2][a1NB2].begin();
                             iB3 !=allDictBondsIdx2D[ha1][ha2][a1NB2].end();
                             iB3++)
                        {
                            for (std::map<ID, std::map<ID, std::map<ID,  
                                 std::vector<aValueSet> > > >::iterator iB4
                                 =iB3->second.begin(); iB4 !=iB3->second.end(); iB4++)
                            {
                                for (std::map<ID, std::map<ID, std::vector<aValueSet> > >::iterator iB5=iB4->second.begin();
                                     iB5 !=iB4->second.end(); iB5++)
                                {
                                    for(std::map<ID, std::vector<aValueSet> >::iterator iB6=iB5->second.begin();
                                          iB6 !=iB5->second.end(); iB6++)
                                    {
                                        for(std::vector<aValueSet>::iterator iB7=iB6->second.begin();
                                            iB7 !=iB6->second.end(); iB7++)
                                        {
                                            tBs2.push_back(*iB7);
                                        }
                                    }
                                }
                            }
                        }
                            
                        aValueSet   tVaS;
                        setValueSet(tVaS, tBs2);
                        
                        if (tVaS.numCodValues > 10 && tVaS.sigValue <0.05)
                        {
                            iB->value        =tVaS.value;
                            iB->valueST      =iB->value;
                            iB->sigValue     =tVaS.sigValue;
                            iB->sigValueST   =iB->sigValue;
                            iB->numCodValues =tVaS.numCodValues;        
                        }
                        else if (ccp4BondsA.find(allAtoms[tPair[0]].ccp4Type) !=ccp4BondsA.end() 
                                 && ccp4BondsA[allAtoms[tPair[0]].ccp4Type].find(allAtoms[tPair[1]].ccp4Type)
                                 !=ccp4BondsA[allAtoms[tPair[0]].ccp4Type].end())
                        {
                            getCCP4Bonds(iB, allAtoms[tPair[0]].ccp4Type, allAtoms[tPair[1]].ccp4Type);
                        }
                        else if (ccp4BondsA.find(allAtoms[tPair[1]].ccp4Type) !=ccp4BondsA.end() 
                                 && ccp4BondsA[allAtoms[tPair[1]].ccp4Type].find(allAtoms[tPair[0]].ccp4Type)
                                 !=ccp4BondsA[allAtoms[tPair[1]].ccp4Type].end())
                        {
                            getCCP4Bonds(iB, allAtoms[tPair[1]].ccp4Type, allAtoms[tPair[0]].ccp4Type);
                                    
                        }
                        else
                        {
                            /*
                            std::vector<BondDict> tBs6;
                            
                            for (std::map<ID, std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                                 std::map<ID, std::map<ID, std::map<ID, int> > > > > > > >::iterator iB2
                                 =allDictBondsIdxD[ha1][ha2][a1NB2].begin();
                                 iB2 !=allDictBondsIdxD[ha1][ha2][a1NB2].end();
                                 iB2++)
                            {
                                for(std::map<ID, std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                                    std::map<ID, std::map<ID, int> > > > > > >::iterator iB3=iB2->second.begin();
                                    iB3!=iB2->second.end(); iB3++)
                                {
                                    for(std::map<ID, std::map<ID, std::map<ID, std::map<ID, std::map<ID,
                                        std::map<ID, int> > > > > >::iterator iB4=iB3->second.begin();
                                        iB4!=iB3->second.end(); iB4++)
                                    {
                                        for(std::map<ID,  std::map<ID, std::map<ID, std::map<ID,
                                            std::map<ID, int> > > > > ::iterator iB5=iB4->second.begin();
                                            iB5!=iB4->second.end(); iB5++)
                                        {
                                            for (std::map<ID,  std::map<ID, std::map<ID,
                                                 std::map<ID, int> > > >::iterator iB6=iB5->second.begin();
                                                 iB6 !=iB5->second.end(); iB6++)
                                            {
                                                for (std::map<ID,  std::map<ID, 
                                                     std::map<ID, int> > >::iterator iB7=iB6->second.begin();
                                                     iB7 != iB6->second.end(); iB7++)
                                                {
                                                    for (std::map<ID,  std::map<ID, int> >::iterator iB8=iB7->second.begin();
                                                            iB8 != iB7->second.end(); iB8++)
                                                    {
                                                        for (std::map<ID, int>::iterator iB9=iB8->second.begin();
                                                                iB9 !=iB8->second.end(); iB9++)
                                                        {
                                                            tBs6.push_back(allDictBonds[iB9->second]);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        
                            dLev = 2;
                            setupTargetBondsUsingSymblDist2(tBs6, iB, as0, as1, dLev);
                            */
                            
                            aValueSet   tVaS;
                            setValueSet(tVaS, tBs2);
                            
                        }
                        
                        std::cout << "iFind 2" << std::endl;
                    }
                }
                else // iFind1 failed a1NB2 
                {
                    /*
                    if (allDictBondsIdx3[ha1][ha2][0].numCodValues > 10 && allDictBondsIdx3[ha1][ha2][0].sigValue <0.04)
                    {
                        iB->value        =allDictBondsIdx3[ha1][ha2][0].value;
                        iB->valueST      =iB->value;
                        iB->sigValue     =allDictBondsIdx3[ha1][ha2][0].sigValue;
                        iB->sigValueST   =iB->sigValue;
                        iB->numCodValues =allDictBondsIdx3[ha1][ha2][0].numCodValues;
                     }
                     */
                                       
                    std::vector<aValueSet> tBs1;
                    
                    for (std::map<ID, std::map<ID, std::map<ID, 
                         std::map<ID, std::map<ID, std::vector<aValueSet> > > > > >::iterator iB1
                         =allDictBondsIdx2D[ha1][ha2].begin();
                         iB1 !=allDictBondsIdx2D[ha1][ha2].end(); iB1++)
                    {
                        for (std::map<ID, std::map<ID, std::map<ID, 
                             std::map<ID, std::vector<aValueSet> > > > >::iterator iB2
                                =iB1->second.begin(); 
                             iB2 !=iB1->second.end(); iB2++)
                        {
                            if (iB2->first == a2NB2)
                            {
                                for (std::map<ID, std::map<ID, std::map<ID,  
                                     std::vector<aValueSet> > > >::iterator iB3
                                     =iB2->second.begin(); iB3 !=iB2->second.end(); iB3++)
                                {
                                    for (std::map<ID, std::map<ID, std::vector<aValueSet> > >::iterator iB4=iB3->second.begin();
                                         iB4 !=iB3->second.end(); iB4++)
                                    {
                                        for(std::map<ID, std::vector<aValueSet> >::iterator iB5=iB4->second.begin();
                                                iB5 !=iB4->second.end(); iB5++)
                                        {
                                            for(std::vector<aValueSet>::iterator iB6=iB5->second.begin();
                                                iB6 !=iB5->second.end(); iB6++)
                                            {
                                                tBs1.push_back(*iB6);
                                                // std::cout << "Bond value " << iB6->value << std::endl;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    aValueSet   tVaS;
                    setValueSet(tVaS, tBs1);
                    
                    // std::cout << "N values " << tVaS.numCodValues << std::endl;
                    // std::cout << "N  siga "  << tVaS.sigValue << std::endl;
                    if (tVaS.numCodValues > 10 && tVaS.sigValue <0.05)
                    {
                        iB->value        =tVaS.value;
                        iB->valueST      =iB->value;
                        iB->sigValue     =tVaS.sigValue;
                        iB->sigValueST   =iB->sigValue;
                        iB->numCodValues =tVaS.numCodValues;
                                
                        // std::cout << "iFind 1 A" << std::endl;           
                    }        
                    else if (ccp4BondsA.find(allAtoms[tPair[0]].ccp4Type) !=ccp4BondsA.end() 
                        && (ccp4BondsA[allAtoms[tPair[0]].ccp4Type].find(allAtoms[tPair[1]].ccp4Type)
                        !=ccp4BondsA[allAtoms[tPair[0]].ccp4Type].end()) )
                    {
                        if(ccp4BondsA[allAtoms[tPair[0]].ccp4Type].find(allAtoms[tPair[1]].ccp4Type)
                            !=ccp4BondsA[allAtoms[tPair[0]].ccp4Type].end())
                        {
                            getCCP4Bonds(iB, allAtoms[tPair[0]].ccp4Type, allAtoms[tPair[1]].ccp4Type);
                        }
                    }
                    else if (ccp4BondsA.find(allAtoms[tPair[1]].ccp4Type) !=ccp4BondsA.end() 
                            && (ccp4BondsA[allAtoms[tPair[1]].ccp4Type].find(allAtoms[tPair[0]].ccp4Type)
                                !=ccp4BondsA[allAtoms[tPair[1]].ccp4Type].end()))
                    {
                        if(ccp4BondsA[allAtoms[tPair[1]].ccp4Type].find(allAtoms[tPair[0]].ccp4Type)
                                !=ccp4BondsA[allAtoms[tPair[0]].ccp4Type].end())
                        {
                            getCCP4Bonds(iB, allAtoms[tPair[1]].ccp4Type, allAtoms[tPair[0]].ccp4Type);
                        }
                    }
                    else
                    {
                        
                        for (std::map<ID, std::map<ID, std::map<ID, 
                         std::map<ID, std::map<ID, std::vector<aValueSet> > > > > >::iterator iB1
                         =allDictBondsIdx2D[ha1][ha2].begin();
                         iB1 !=allDictBondsIdx2D[ha1][ha2].end(); iB1++)
                        {
                            for (std::map<ID, std::map<ID, std::map<ID, 
                                 std::map<ID, std::vector<aValueSet> > > > >::iterator iB2
                                   =iB1->second.begin(); 
                                 iB2 !=iB1->second.end(); iB2++)
                            {
                                if (iB2->first != a2NB2)
                                {
                                    for (std::map<ID, std::map<ID, std::map<ID,  
                                         std::vector<aValueSet> > > >::iterator iB3
                                         =iB2->second.begin(); iB3 !=iB2->second.end(); iB3++)
                                    {
                                        for (std::map<ID, std::map<ID, std::vector<aValueSet> > >::iterator iB4=iB3->second.begin();
                                             iB4 !=iB3->second.end(); iB4++)
                                        {
                                            for(std::map<ID, std::vector<aValueSet> >::iterator iB5=iB4->second.begin();
                                                iB5 !=iB4->second.end(); iB5++)
                                            {
                                                for(std::vector<aValueSet>::iterator iB6=iB5->second.begin();
                                                    iB6 !=iB5->second.end(); iB6++)
                                                {
                                                    tBs1.push_back(*iB6);
                                                    // std::cout << "Bond value " << iB6->value << std::endl;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    
                        aValueSet   tVaS1;
                        setValueSet(tVaS1, tBs1);
                        /*
                        dLev = 1;
                        setupTargetBondsUsingSymblDist2(tBs6, iB, as0, as1, dLev);
                         */
                        iB->value        =tVaS1.value;
                        iB->valueST      =iB->value;
                        iB->sigValue     =tVaS1.sigValue;
                        iB->sigValueST   =iB->sigValue;
                        iB->numCodValues =tVaS1.numCodValues;
                        
                        std::cout << "iFind 1" << std::endl;
                    }
                }
                
            }
            else if (ccp4BondsA.find(allAtoms[tPair[0]].ccp4Type) !=ccp4BondsA.end() 
                            && (ccp4BondsA[allAtoms[tPair[0]].ccp4Type].find(allAtoms[tPair[1]].ccp4Type)
                                !=ccp4BondsA[allAtoms[tPair[0]].ccp4Type].end()
                               || ccp4BondsA[allAtoms[tPair[0]].ccp4Type].find(".") !=ccp4BondsA[allAtoms[tPair[0]].ccp4Type].end()))
            {
                // Without hash code matching for those combination of two atoms. try CCP4 type matching
                // It is very unlikely we have ccp4 type matching when hash code matching failed.
                
                if(ccp4BondsA[allAtoms[tPair[0]].ccp4Type].find(allAtoms[tPair[1]].ccp4Type)
                                !=ccp4BondsA[allAtoms[tPair[0]].ccp4Type].end())
                {
                    
                    getCCP4Bonds(iB, allAtoms[tPair[0]].ccp4Type, allAtoms[tPair[1]].ccp4Type);
                }
                else
                {
                    
                    
                    getCCP4Bonds(iB, allAtoms[tPair[0]].ccp4Type, ".");
                }
            }
            else if (ccp4BondsA.find(allAtoms[tPair[1]].ccp4Type) !=ccp4BondsA.end() 
                            && (ccp4BondsA[allAtoms[tPair[1]].ccp4Type].find(allAtoms[tPair[0]].ccp4Type)
                                !=ccp4BondsA[allAtoms[tPair[1]].ccp4Type].end()
                                || ccp4BondsA[allAtoms[tPair[1]].ccp4Type].find(".") !=ccp4BondsA[allAtoms[tPair[1]].ccp4Type].end()))
            {
                // Without hash code matching for those combination of two atoms. try CCP4 type matching
                // It is very unlikely we have ccp4 type matching when hash code matching failed.
                if(ccp4BondsA[allAtoms[tPair[1]].ccp4Type].find(allAtoms[tPair[0]].ccp4Type)
                                !=ccp4BondsA[allAtoms[tPair[1]].ccp4Type].end())
                {
                    getCCP4Bonds(iB, allAtoms[tPair[1]].ccp4Type, allAtoms[tPair[0]].ccp4Type);
                }
                else
                {
                    getCCP4Bonds(iB, allAtoms[tPair[1]].ccp4Type, ".");
                }
            }
            else
            {
                std::cout << "Could not find the basic bond classes in Cod and CCP4-energy-type" << std::endl
                          << "for the dictionary bond of atoms " <<iB->atoms[0] << " and "
                          << iB->atoms[1] <<std::endl;
                
                for (std::map<ID, int>::iterator iAt=iB->fullAtoms.begin();
                        iAt != iB->fullAtoms.end(); iAt++)
                {
                    std::cout << "chemType : " << allAtoms[iAt->second].chemType << " atom COD classes : " 
                              << allAtoms[iAt->second].codClass << std::endl;
                }
            }
            
            
            //std::cout << "The bond value is now " << iB->value << std::endl;
            //if (iB->hasCodValue)
            //{
            //    std::cout << "It have exact matches to COD atoms " << std::endl;
            //}
            //else
            //{
            //    std::cout << "The bond value is generated by " << std::endl
            //              << "either averaging or shortest distance to COD"
            //              << std::endl;
            //}
    }
    
    void CodClassify::getCCP4Bonds(std::vector<BondDict>::iterator tB, 
                                   ID tAtom1, ID tAtom2)
    {
        
        
        std::string aOrdS = tB->order.substr(0,4);
        StrUpper(aOrdS);
        
        std::string tAtm1, tAtm2; 
        
        if (tAtom1.find(".") == std::string::npos)
        {
            tAtm1 = tAtom1;
            tAtm2 = tAtom2;
        }
        else
        {
            tAtm1 = tAtom2;
            tAtm2 = tAtom1;
        }
        
        if (ccp4BondsA.find(tAtm1) != ccp4BondsA.end() )
        {
            if (ccp4BondsA[tAtm1].find(tAtm2) != ccp4BondsA[tAtm1].end())
            {
                
                        
                if (ccp4BondsA[tAtm1][tAtm2].find(aOrdS) !=ccp4BondsA[tAtm1][tAtm2].end())
                { 
                   
                    tB->value    = ccp4BondsA[tAtm1][tAtm2][aOrdS]["length"];
                    tB->sigValue = ccp4BondsA[tAtm1][tAtm2][aOrdS]["sigValue"];
                }
                else
                {
                    tB->value    = 
                    ccp4BondsA[tAtm1][tAtm2][ccp4BondsA[tAtm1][tAtm2].begin()->first]["length"];
                    tB->sigValue = 
                    ccp4BondsA[tAtm1][tAtm2][ccp4BondsA[tAtm1][tAtm2].begin()->first]["sigValue"];
                }
            }
            else if (ccp4BondsA[tAtm1].find(".") != ccp4BondsA[tAtm1].end())
            {
                if (ccp4BondsA[tAtm1]["."].find(aOrdS) !=ccp4BondsA[tAtm1]["."].end())
                { 
                    tB->value    = ccp4BondsA[tAtm1]["."][aOrdS]["length"];
                    tB->sigValue = ccp4BondsA[tAtm1]["."][aOrdS]["sigValue"];
                }
                else
                {
                    tB->value    = 
                    ccp4BondsA[tAtm1]["."][ccp4BondsA[tAtm1]["."].begin()->first]["length"];
                    tB->sigValue = 
                    ccp4BondsA[tAtm1]["."][ccp4BondsA[tAtm1]["."].begin()->first]["sigValue"];
                }
            }
            else
            {
                std::cout << "Bug: can not find the bond type beginning with CCP4 atom type " 
                          << tAtm1  << std::endl;
                exit(1);
            }
        }
        else if (tAtm2.compare(".") != 0 && ccp4BondsA.find(tAtm2) != ccp4BondsA.end() )
        {
            if (ccp4BondsA[tAtm2].find(tAtm1) != ccp4BondsA[tAtm2].end())
            {
                if (ccp4BondsA[tAtm2][tAtm1].find(aOrdS) !=ccp4BondsA[tAtm2][tAtm1].end())
                { 
                    tB->value    = ccp4BondsA[tAtm2][tAtm1][aOrdS]["length"];
                    tB->sigValue = ccp4BondsA[tAtm2][tAtm1][aOrdS]["sigValue"];
                }
                else
                {
                    tB->value    = 
                    ccp4BondsA[tAtm2][tAtm1][ccp4BondsA[tAtm2][tAtm1].begin()->first]["length"];
                    tB->sigValue = 
                    ccp4BondsA[tAtm2][tAtm1][ccp4BondsA[tAtm2][tAtm1].begin()->first]["sigValue"];
                }
            }
            else if (ccp4BondsA[tAtm2].find(".") != ccp4BondsA[tAtm2].end())
            {
                if (ccp4BondsA[tAtm2]["."].find(aOrdS) !=ccp4BondsA[tAtm2]["."].end())
                { 
                    tB->value    = ccp4BondsA[tAtm2]["."][aOrdS]["length"];
                    tB->sigValue = ccp4BondsA[tAtm2]["."][aOrdS]["sigValue"];
                }
                else
                {
                    tB->value    = 
                    ccp4BondsA[tAtm2]["."][ccp4BondsA[tAtm2]["."].begin()->first]["length"];
                    tB->sigValue = 
                    ccp4BondsA[tAtm2]["."][ccp4BondsA[tAtm2]["."].begin()->first]["sigValue"];
                }
            }
            else
            {
                std::cout << "Bug: can not find the bond type beginning with CCP4 atom type " 
                          << tAtm2  << std::endl;
                exit(1);
            }
        }
    }
    
    
    void CodClassify::getCCP4Bonds2(std::vector<BondDict>::iterator tB, 
                                    ID tAtom1, ID tAtom2)
    {
        
        std::string aOrdS = tB->order.substr(0,4);
        StrUpper(aOrdS);
        
        std::string tAtm1, tAtm2; 
        
        if (tAtom1.find(".") == std::string::npos)
        {
            tAtm1 = tAtom1;
            tAtm2 = tAtom2;
        }
        else
        {
            tAtm1 = tAtom2;
            tAtm2 = tAtom1;
        }
        
        std::map<int, bool> posM;
        posM[111]=false;
        posM[110]=false;
        posM[100]=false;
        posM[121]=false;
        posM[120]=false;
        posM[211]=false;
        posM[210]=false;
        posM[200]=false;
        posM[221]=false;
        posM[220]=false;
        
        posM[000]=false;
        
        // Different situations
        
        if (ccp4BondsA.find(tAtm1) != ccp4BondsA.end())
        {
            if (ccp4BondsA[tAtm1].find(tAtm2) != ccp4BondsA[tAtm1].end())
            {
                if(ccp4BondsA[tAtm1][tAtm2].find(aOrdS) !=ccp4BondsA[tAtm1][tAtm2].end())
                {
                    // type1, type2 and order are matched
                    posM[111]=true;
                }
                posM[110] = true;
            }
            else if (ccp4BondsA[tAtm1].find(".") != ccp4BondsA[tAtm1].end())
            {
                if (ccp4BondsA[tAtm1]["."].find(aOrdS) !=ccp4BondsA[tAtm1]["."].end())
                {
                    posM[121] = true;
                }
                else 
                {
                    posM[120] = true;
                }
            }
            posM[100] = true;
        }
        else if (ccp4BondsA.find(tAtm2) != ccp4BondsA.end())
        {
            if (ccp4BondsA[tAtm2].find(tAtm1) != ccp4BondsA[tAtm2].end())
            {
                if(ccp4BondsA[tAtm2][tAtm1].find(aOrdS) !=ccp4BondsA[tAtm2][tAtm1].end())
                {
                    // type1, type2 and order are matched
                    posM[211]=true;
                }
                posM[210] = true;
            }
            else if (ccp4BondsA[tAtm2].find(".") != ccp4BondsA[tAtm2].end())
            {
                if (ccp4BondsA[tAtm2]["."].find(aOrdS) !=ccp4BondsA[tAtm2]["."].end())
                {
                    posM[221] = true;
                }
                else 
                {
                    posM[220] = true;
                }
            }
            posM[200] = true;
        }
        else
        {
            posM[000] = true;
        }
        
        setCCP4BondByMode(tB, tAtm1, tAtm2, aOrdS, posM);
        
    }

    void CodClassify::setCCP4BondByMode(std::vector<BondDict>::iterator tB,
                                        ID tAtm1, ID tAtm2, ID tOrdS,
                                        std::map<int, bool> & tPosM)
    {
        if (tPosM[111])
        {
            tB->value    = ccp4BondsA[tAtm1][tAtm2][tOrdS]["length"];
            tB->sigValue = ccp4BondsA[tAtm1][tAtm2][tOrdS]["sigValue"];
            tB->nLevel   = 0;
        }
        else if (tPosM[211])
        {
            tB->value    = ccp4BondsA[tAtm2][tAtm1][tOrdS]["length"];
            tB->sigValue = ccp4BondsA[tAtm2][tAtm1][tOrdS]["sigValue"];
            tB->nLevel   = 0;
        }
        else if (tPosM[110])
        {
            setCCP4Mix(tB, tAtm1, tAtm2, tOrdS, 110);
        }
        else if(tPosM[210])
        {
            setCCP4Mix(tB, tAtm1, tAtm2, tOrdS, 210);
        }
        else if(tPosM[121] && tPosM[221])
        {
            setCCP4Mix(tB, tAtm1, tAtm2, tOrdS, 342);
        }
        else if(tPosM[121] && tPosM[220])
        {
            setCCP4Mix(tB, tAtm1, tAtm2, tOrdS, 341);
        }
    }
    
    void CodClassify::setCCP4Mix(std::vector<BondDict>::iterator tB,
                               ID tAtm1, ID tAtm2, 
                               ID tOrdS, int tMode)
    {
    }
    
    void CodClassify::groupCodMetBonds()
    {
        try
        {
            // std::cout << "Clustering all COD Metal bonds " << std::endl;
            
            // should be something like std::string tNewCodBondFileName(clibMonDir + "/list/bonds.txt");
            //std::string clibMonDir(std::getenv("CLIBD_MON"));
            //std::string fName = clibMonDir+"allMetalBonds.table";
            // std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
            std::string fName = libmolTabDir  + "/allMetalBonds.table";
            std::ifstream codMetBondFile(fName.c_str());
            if(codMetBondFile.is_open())
            {
                
                
                std::string tRecord="";
                
                while(!codMetBondFile.eof())
                {
                    std::getline(codMetBondFile, tRecord);
                    tRecord = TrimSpaces(tRecord);
                    std::vector<std::string> tBuf;
                    StrTokenize(tRecord, tBuf);
                    
                    for (int i=0; i < (int)tBuf.size(); i++)
                    {
                        tBuf[i] = TrimSpaces(tBuf[i]);
                    }
                    
                    if ((int)tBuf.size() ==11)
                    {
                           
                        int ha2 = StrToInt(tBuf[1]);
                        int ha3 = StrToInt(tBuf[3]);
                        
                        //std::cout << tBuf[0] << " and " << tBuf[2] << " and " 
                        //         << tBuf[7] << std::endl;
                        //std::cout << "ha2 " << ha2 << " ha3 " << ha3 << std::endl;
                        allDictPreMetBonds[tBuf[0]][tBuf[2]][ha2][ha3]["value"]     
                                = StrToReal(tBuf[4]);
                        allDictPreMetBonds[tBuf[0]][tBuf[2]][ha2][ha3]["sigma"]     
                                = StrToReal(tBuf[5]);
                        allDictPreMetBonds[tBuf[0]][tBuf[2]][ha2][ha3]["numValues"] 
                                = StrToReal(tBuf[6]);
                        
                        
                        allDictMetBonds[tBuf[0]][tBuf[2]][ha2][ha3][tBuf[7]]["value"]     
                                = StrToReal(tBuf[8]);
                        allDictMetBonds[tBuf[0]][tBuf[2]][ha2][ha3][tBuf[7]]["sigma"]     
                                = StrToReal(tBuf[9]);
                        allDictMetBonds[tBuf[0]][tBuf[2]][ha2][ha3][tBuf[7]]["numValues"]     
                                = StrToReal(tBuf[10]);
                   
                    }
                    else if ((int)tBuf.size() ==10)
                    {
                        int ha2 = StrToInt(tBuf[1]);
                        int ha3 = StrToInt(tBuf[3]);
                        std::string tS= "NONE";
                        //std::cout << tBuf[0] << " and " << tBuf[2] << std::endl;
                        //std::cout << "ha2 " << ha2 << " ha3 " << ha3 << std::endl;
                        allDictPreMetBonds[tBuf[0]][tBuf[2]][ha2][ha3]["value"]     
                                = StrToReal(tBuf[4]);
                        allDictPreMetBonds[tBuf[0]][tBuf[2]][ha2][ha3]["sigma"]     
                                = StrToReal(tBuf[5]);
                        allDictPreMetBonds[tBuf[0]][tBuf[2]][ha2][ha3]["numValues"] 
                                = StrToReal(tBuf[6]);
                        
                        
                        allDictMetBonds[tBuf[0]][tBuf[2]][ha2][ha3][tS]["value"]     
                                = StrToReal(tBuf[7]);
                        allDictMetBonds[tBuf[0]][tBuf[2]][ha2][ha3][tS]["sigma"]     
                                = StrToReal(tBuf[8]);
                        allDictMetBonds[tBuf[0]][tBuf[2]][ha2][ha3][tS]["numValues"]     
                                = StrToReal(tBuf[9]);
                    }
                }
                codMetBondFile.close();
            }
        }
        catch (std::exception & e)
        {
            std::cout << e.what() << std::endl;
        }
        std::cout << "Group COD bonds involving metal elements " << std::endl;
    }
    
    
    
    void CodClassify::searchCodMetBonds(std::vector<BondDict>::iterator iMB)
    {
        std::vector<int> tPair;
        for (std::map<ID, int>::iterator iA=iMB->fullAtoms.begin();
                    iA !=iMB->fullAtoms.end(); iA++)
        {
            tPair.push_back(iA->second);
        }
            
        
        int ha1, ha2, tIdx=-1;
        ID  id1, id2,  a2NB;
        std::vector<int> nonMetNB2;
        
        if(allAtoms[tPair[0]].isMetal) 
        {
            tIdx = tPair[1];
            id1   =allAtoms[tPair[0]].chemType;
            id2   =allAtoms[tPair[1]].chemType;
            
            
            ha1   = (int)allAtoms[tPair[0]].connAtoms.size();
            for (std::vector<int>::iterator iN=allAtoms[tPair[1]].connAtoms.begin();
                     iN != allAtoms[tPair[1]].connAtoms.end(); iN++)
            {
                if (!allAtoms[*iN].isMetal)
                {
                    nonMetNB2.push_back(*iN);
                }
            }
            ha2   = (int)nonMetNB2.size();
           
        }
        else if (allAtoms[tPair[1]].isMetal)
        {
            tIdx = tPair[0];
            id1   =allAtoms[tPair[1]].chemType;
            id2   =allAtoms[tPair[0]].chemType;
            
            ha1   = (int)allAtoms[tPair[1]].connAtoms.size();
            for (std::vector<int>::iterator iN=allAtoms[tPair[0]].connAtoms.begin();
                     iN != allAtoms[tPair[0]].connAtoms.end(); iN++)
            {
                if (!allAtoms[*iN].isMetal)
                {
                    nonMetNB2.push_back(*iN);
                }
            }
            ha2   = (int)nonMetNB2.size();
        }
        else
        {
            std::cout << "a bond contains no metal atom " << std::endl;
            exit(1);
        }
        
        
        if (tIdx >=0)
        {
            std::list<ID> tIds;
            for (std::vector<int>::iterator iN= nonMetNB2.begin();
                 iN != nonMetNB2.end(); iN++)
            {
                tIds.push_back(allAtoms[*iN].chemType);
            }
            tIds.sort(compareNoCase);
            
            for (std::list<ID>::iterator iID = tIds.begin(); iID !=tIds.end();
                  iID++)
            {
                a2NB.append("_"+ *iID);
            }
        }
        
        
        std::cout << "element : " <<id1 << "  CN  " << ha1 <<  std::endl
                  << "element : " <<id2 << "  CN  " << ha2 << " NB " 
                  << a2NB <<  std::endl;
          
        std::map<ID, std::map<ID, std::map<int, std::map<int, 
                 std::map<ID, std::map<ID, REAL> > > > > >::iterator iFind1 = allDictMetBonds.find(id1);
        if (iFind1 !=allDictMetBonds.end())
        {
            std::map<ID, std::map<int, std::map<int, 
            std::map<ID, std::map<ID, REAL> > > > >::iterator iFind2 = allDictMetBonds[id1].find(id2);
            if (iFind2 != allDictMetBonds[id1].end())
            {
                std::map<int, std::map<int, 
                std::map<ID, std::map<ID, REAL> > > >::iterator iFind3 = allDictMetBonds[id1][id2].find(ha1);
                if (iFind3 != allDictMetBonds[id1][id2].end())
                {
                    std::map<int, std::map<ID, std::map<ID, REAL> > >::iterator iFind4
                    = allDictMetBonds[id1][id2][ha1].find(ha2);
                    if (iFind4 != allDictMetBonds[id1][id2][ha1].end())
                    {
                        std::map<ID, std::map<ID, REAL> >::iterator iFind5 
                                = allDictMetBonds[id1][id2][ha1][ha2].find(a2NB);
                        if (iFind5 != allDictMetBonds[id1][id2][ha1][ha2].end())
                        {
                            if (allDictMetBonds[id1][id2][ha1][ha2][a2NB]["numValues"] > 5.0)
                            {
                                iMB->value    = allDictMetBonds[id1][id2][ha1][ha2][a2NB]["value"];
                                iMB->sigValue = allDictMetBonds[id1][id2][ha1][ha2][a2NB]["sigma"];
                                iMB->numCodValues = (int)allDictMetBonds[id1][id2][ha1][ha2][a2NB]["numValues"];
                            }
                            else
                            {
                                // Using the parent bond value set
                                iMB->value    = allDictPreMetBonds[id1][id2][ha1][ha2]["value"];
                                iMB->sigValue = allDictPreMetBonds[id1][id2][ha1][ha2]["sigma"];
                                iMB->numCodValues = (int)allDictPreMetBonds[id1][id2][ha1][ha2]["numValues"];
                            }
                        }
                        else
                        {
                            // Using the parent bond value set
                            iMB->value    = allDictPreMetBonds[id1][id2][ha1][ha2]["value"];
                            iMB->sigValue = allDictPreMetBonds[id1][id2][ha1][ha2]["sigma"];
                            iMB->numCodValues = allDictPreMetBonds[id1][id2][ha1][ha2]["numValues"];
                        }
                    }
                    else
                    {
                        std::cout << "COD does not have bond between ( " << id1 
                                  << " coordination number " << ha1 << " and "
                                  << id2 << " coordination number " << ha2 
                                  << " )" << std::endl;
                        exit(1);
                    }
                }
                else
                {
                    ID  aID2="", aID2NB="";
                    int aIdx2=-1, tDIdx=12, dIdx=12;
                    std::map<ID, std::map<int, std::map<int, 
                    std::map<ID, std::map<ID, REAL> > > > >::iterator iId2;
                    for (iId2 =allDictMetBonds[id1].begin(); 
                         iId2 !=allDictMetBonds[id1].end(); iId2++)
                    {
                        std::map<int, std::map<int, 
                        std::map<ID, std::map<ID, REAL> > > >::iterator iCN1;
                        for (iCN1=iId2->second.begin(); iCN1 !=iId2->second.end(); 
                             iCN1++)
                        {
                            std::map<int, std::map<ID, 
                            std::map<ID, REAL> > >::iterator iCN2;
                            for (iCN2 = iCN1->second.begin(); iCN2 !=iCN1->second.end();
                                    iCN2++)
                            {
                                if (iCN1->first==ha1)
                                {
                                    tDIdx = abs(iCN2->first-ha2);
                                    if(tDIdx < dIdx)
                                    {
                                        aID2  = iId2->first;
                                        aIdx2 = iCN2->first;
                                    }
                                }  
                            }   
                        }
                    }
                    if (aIdx2 !=-1)
                    {
                        std::vector<std::map<ID,REAL> > tBondVect;
                        for (std::map<ID, std::map<ID, REAL> >::iterator iId2NB=
                                allDictMetBonds[id1][aID2][ha1][aIdx2].begin();
                                iId2NB != allDictMetBonds[id1][aID2][ha1][aIdx2].end();
                                iId2NB++)
                        {
                            // Using the parent bond value set
                            tBondVect.push_back(iId2NB->second);
                            setupTargetMetBondsUsingMean(tBondVect, iMB);
                        }
                    }
                    else
                    {
                        // searchSimMetBonds(id1Keys, iMB);
                        std::cout << id1 << " does not have coordination number " << ha1 << std::endl;
                        exit(1);
                    }
                }
            }
            else
            {
                std::cout << " No COD bond data between element " << id1 
                        << " and element " << id2 << std::endl;
                exit(1);
            }
            
        }
        else
        {
            std::cout << "Could not find metal element " << id1 << std::endl;
            exit(1);
        }
           
    }
    
    void CodClassify::setValueSet(aValueSet& tVs, std::vector<aValueSet>& tVecVs)
    {
        REAL aSum =0.0, sum1 =0.0, sum2=0.0;
        int  nSum = 0;
        
        tVs.numCodValues =0;
        tVs.value        =0.0;
        tVs.sigValue     =0.0;
        
        for (std::vector<aValueSet>::iterator iVs=tVecVs.begin();
                iVs != tVecVs.end(); iVs++)
        {
            aSum +=(iVs->value*iVs->numCodValues);
            sum1 +=((iVs->numCodValues-1)*iVs->sigValue*iVs->sigValue
                     + iVs->numCodValues*iVs->value*iVs->value);
            nSum += iVs->numCodValues;
        }
        
        if (nSum !=0)
        {
            tVs.value = aSum/nSum;
            tVs.numCodValues = nSum;
            sum2 = tVs.value*aSum;
            
        }
        else
        {
            std::cout << "Bug: Acedrg database contains zero observation term "
                      << std::endl;
            exit(1);
        }
        
        if (nSum > 1)
        {
            tVs.sigValue = sqrt(fabs(sum1-sum2)/(nSum-1));
        }
        else if (nSum==1)
        {
            tVs.sigValue = sqrt(fabs(sum1-sum2)/nSum);
        }
        else
        {
            std::cout << "Bug: Acedrg database contains zero observation term "
                      << std::endl;
            exit(1);
        }
    }
    
    
    void CodClassify::setupTargetMetBondsUsingMean(std::vector<std::map<ID,REAL> > tBondVect,  
                                        std::vector<BondDict>::iterator iMB)
    {
        REAL tSumValues=0.0, tSumSqValues=0.0, tMeanSqValue=0.0;
        REAL tMeanValue=0.0, tSqMeanValue=0.0, tStVar=0.0;
        int nTot=0;
        for (int i=0; i < (int)tBondVect.size(); i++)
        {
            nTot +=((int)tBondVect[i]["numValues"]);
            tSumValues +=(tBondVect[i]["value"]*tBondVect[i]["numValues"]);
            tSumSqValues+=(tBondVect[i]["value"]*tBondVect[i]["value"]*tBondVect[i]["numValues"]);
        }
        
        tMeanValue    = tSumValues/nTot;
        
        tSqMeanValue = (tMeanValue*tMeanValue);
        
        tMeanSqValue  = (tSumSqValues/nTot);
        
        tStVar  = std::fabs(tMeanSqValue-tSqMeanValue);
        
        tStVar  = std::sqrt(tStVar);
     
        iMB->value        = tMeanValue;
        iMB->sigValue     = tStVar;
        
        if(iMB->sigValue < 0.02)
        {
            iMB->sigValue = 0.02;
        }
        
        iMB->numCodValues = nTot;
    }
    
    void CodClassify::setupTargetBondsUsingSymblDist(std::vector<BondDict>& tBonds,
                                                     std::vector<BondDict>::iterator tB,
                                                     int tLev)
    {
        
        std::cout << "Using the matched value by the distance " << std::endl;
        // search the closest cod-classes
        int shortestDist =1000000, iPos=-1;
        
        int dist1   =1000000;
        int dist2   =1000000;
        
        
        std::vector<std::vector<ID> > targetNBs;
        
        for (std::map<ID,int>::iterator iBT=tB->fullAtoms.begin();
                iBT !=tB->fullAtoms.end(); iBT++)
        {
            std::vector<ID> tV;
            
            StrTokenize(allAtoms[iBT->second].codNBSymb, tV, ':');
            targetNBs.push_back(tV);
        }
        
        for (int ib=0; ib < (int)tBonds.size(); ib++)
        {
            std::vector<ID> tIdV0, tIdV1;
            
            StrTokenize(tBonds[ib].atomsNBRep[0], tIdV0, ':');
            StrTokenize(tBonds[ib].atomsNBRep[1], tIdV1, ':');
            //std::cout << tBonds[ib].atomsNBRep[0] << std::endl
            //      << tBonds[ib].atomsNBRep[1] << std::endl;
            
            if ((int)tIdV0.size() == (int)targetNBs[0].size() &&
                 (int)tIdV1.size() == (int)targetNBs[1].size()   )
            {
                
                dist1 = codAtomsDist(targetNBs[0], tIdV0, tLev) 
                      + codAtomsDist(targetNBs[1], tIdV1, tLev);
            }
            if ((int)tIdV0.size() == (int)targetNBs[1].size() &&
                 (int)tIdV1.size() == (int)targetNBs[0].size()   )
            {
                
                dist2 = codAtomsDist(targetNBs[0], tIdV1, tLev) 
                      + codAtomsDist(targetNBs[1], tIdV0, tLev);
            }                  
            
            if (dist1 < shortestDist )
            {
                shortestDist=dist1;
                iPos = ib;
            }
            if (dist2 < shortestDist )
            {
                shortestDist=dist2;
                iPos = ib;
            }
        }
        
        if (iPos !=-1)
        {
           tB->value        = tBonds[iPos].value;
           tB->sigValue     = tBonds[iPos].sigValue;
           if (tB->sigValue < 0.02)
           {
               tB->sigValue = 0.02;
           }
           tB->numCodValues = tBonds[iPos].numCodValues; 
           
           //std::cout << "target bond atom cod-classes: " << std::endl;
           //for (std::map<ID,int>::iterator iA = tB->fullAtoms.begin();
           //     iA != tB->fullAtoms.end(); iA++)
           //{
           //   std::cout << allAtoms[iA->second].codClass << std::endl;
           //}
        
           //std::cout << "found the closest bond atom cod-classes: " << std::endl
           //          << tBonds[iPos].atomsCodClasses[0] << std::endl 
           //          << tBonds[iPos].atomsCodClasses[1] << std::endl;
        }
        else
        {
           //std::cout << "target bond atom cod-classes: " << std::endl;
           //for (std::map<ID,int>::iterator iA = tB->fullAtoms.begin();
           //     iA != tB->fullAtoms.end(); iA++)
           //{
            //  std::cout << allAtoms[iA->second].codClass << std::endl;
            //  std::cout << allAtoms[iA->second].hashingValue << std::endl;
          // }
           
           //std::cout << iPos << std::endl;
           //std::cout << "Search bond atom cod-classes: " << std::endl;
           
           for (int ib1=0; ib1 < (int)tBonds.size(); ib1++)
           {
               std::cout << "Bond " << ib1 << std::endl;
               std::cout << "atom 1: " << std::endl 
                         << "COD code " << tBonds[ib1].atomsCodClasses[0] << std::endl
                         << "atom 2: " << tBonds[ib1].atomsCodClasses[1] << std::endl;
           }
           exit(1);
        }       
    }
    
    void CodClassify::setupTargetBondsUsingSymblDist2(std::vector<BondDict>& tBonds,
                                                     std::vector<BondDict>::iterator tB,
                                                     int tAs0, int tAs1, int tLev)
    {
        int a0Idx=0, a1Idx=0;
        int iPos = -1, mDiff=1000000;
        
        if (tLev==2)
        {
            // Last two NB2 
            // the last codNB2Symb
            std::vector<ID> tV1;
            StrTokenize(allAtoms[tB->atomsIdx[tAs1]].codNB2Symb, tV1, ':'); 
            for (int i=0; i < (int)tV1.size(); i++)
            {
                a1Idx+=StrToInt(tV1[i]);
            }
            
             // scan the group of angles
            for (int i=0; i <(int)tBonds.size(); i++)
            {
                std::vector<ID> tAV1;
                int tA1Idx=0; 
                StrTokenize(tBonds[i].atomsNB2Rep[1], tAV1, ':');
                for (int j=0; j < (int)tV1.size(); j++)
                {
                    tA1Idx+=StrToInt(tV1[j]);
                }
                
                int tDiff= abs(tA1Idx-a1Idx);
                
                if (tDiff < mDiff)
                {
                    iPos = i;
                    mDiff = tDiff;
                }
            }            
        }
        if (tLev==1)
        {
            // Last two NB2 
            // the last codNB2Symb
            std::vector<ID> tV0, tV1;
            StrTokenize(allAtoms[tB->atomsIdx[tAs0]].codNB2Symb, tV0, ':'); 
            for (int i=0; i < (int)tV0.size(); i++)
            {
                a0Idx+=StrToInt(tV0[i]);
            }
            StrTokenize(allAtoms[tB->atomsIdx[tAs1]].codNB2Symb, tV1, ':'); 
            for (int i=0; i < (int)tV1.size(); i++)
            {
                a1Idx+=StrToInt(tV1[i]);
            }
            
             // scan the group of angles
            for (int i=0; i <(int)tBonds.size(); i++)
            {
                std::vector<ID> tAV0;
                int tA0Idx=0; 
                StrTokenize(tBonds[i].atomsNB2Rep[0], tAV0, ':');
                for (int j=0; j < (int)tV0.size(); j++)
                {
                    tA0Idx+=StrToInt(tV0[j]);
                }
                
                std::vector<ID> tAV1;
                int tA1Idx=0; 
                StrTokenize(tBonds[i].atomsNB2Rep[1], tAV1, ':');
                for (int j=0; j < (int)tV1.size(); j++)
                {
                    tA1Idx+=StrToInt(tV1[j]);
                }
                
                int tDiff= abs(tA0Idx-a0Idx)+abs(tA1Idx-a1Idx);
                
                if (tDiff < mDiff)
                {
                    iPos = i;
                    mDiff = tDiff;
                }
            }            
        }
        
        if (iPos !=-1)
        {
           tB->value        = tBonds[iPos].value;
           tB->sigValue     = tBonds[iPos].sigValue;
           if (tB->sigValue < 0.02)
           {
               tB->sigValue = 0.02;
           }
           tB->numCodValues = tBonds[iPos].numCodValues; 
           
           //std::cout << "target bond atom cod-classes: " << std::endl;
           //for (std::map<ID,int>::iterator iA = tB->fullAtoms.begin();
           //     iA != tB->fullAtoms.end(); iA++)
           //{
           //   std::cout << allAtoms[iA->second].codClass << std::endl;
           //}
        
           //std::cout << "found the closest bond atom cod-classes: " << std::endl
           //          << tBonds[iPos].atomsCodClasses[0] << std::endl 
           //          << tBonds[iPos].atomsCodClasses[1] << std::endl;
        }
        else
        {
           //std::cout << "target bond atom cod-classes: " << std::endl;
           //for (std::map<ID,int>::iterator iA = tB->fullAtoms.begin();
           //     iA != tB->fullAtoms.end(); iA++)
           //{
            //  std::cout << allAtoms[iA->second].codClass << std::endl;
            //  std::cout << allAtoms[iA->second].hashingValue << std::endl;
          // }
           
           //std::cout << iPos << std::endl;
           //std::cout << "Search bond atom cod-classes: " << std::endl;
           
           for (int ib1=0; ib1 < (int)tBonds.size(); ib1++)
           {
               std::cout << "Bond " << ib1 << std::endl;
               std::cout << "atom 1: " << std::endl 
                         << "COD code " << tBonds[ib1].atomsCodClasses[0] << std::endl
                         << "atom 2: " << tBonds[ib1].atomsCodClasses[1] << std::endl;
           }
           exit(1);
        }     
    }
    void CodClassify:: setupTargetBondsUsingMean(std::vector<BondDict> & tBonds,
                                   std::vector<BondDict>::iterator tB)
    {
        std::cout << "Using the meaning value " << std::endl;
        REAL tSumValues=0.0, tSumSqValues=0.0, tMeanSqValue=0.0;
        REAL tMeanValue=0.0, tSqMeanValue=0.0, tStVar=0.0;
        int nTot=0;
        for (std::vector<BondDict>::iterator iB = tBonds.begin(); 
                iB !=tBonds.end(); iB++)
        {
            nTot +=iB->numCodValues;
            tSumValues +=(iB->value*(iB->numCodValues));
            tSumSqValues+=(iB->value*(iB->value)*(iB->numCodValues));
        }
        
        tMeanValue    = tSumValues/nTot;
        
        tSqMeanValue = (tMeanValue*tMeanValue);
        
        tMeanSqValue  = (tSumSqValues/nTot);
        
        tStVar  = std::fabs(tMeanSqValue-tSqMeanValue);
        
        tStVar  = std::sqrt(tStVar);
     
        tB->value        = tMeanValue;
        tB->sigValue     = tStVar;
        
        if(tB->sigValue < 0.01)
        {
            tB->sigValue = 0.01;
        }
        
        tB->numCodValues = nTot;
    }
    
    void CodClassify::setupTargetBondsUsingValueSetMean(std::vector<aValueSet>& tSets, 
                                                        std::vector<BondDict>::iterator tB)
    {
        if (tSets.size() !=0)
        {
            REAL aSum=0.0, vSum1=0.0, vSum2=0.0, vSum3=0.0;
            REAL aCount = 0;
            for (std::vector<aValueSet>::iterator iSet=tSets.begin();
                    iSet != tSets.end(); iSet++)
            {
                aSum+=(iSet->value*iSet->numCodValues);
                aCount+=iSet->numCodValues;
                vSum1 +=((iSet->numCodValues-1)*iSet->sigValue*iSet->sigValue 
                          + iSet->numCodValues*iSet->value*iSet->value);
                vSum2 +=(iSet->value*iSet->numCodValues);
            }
            
            if (aCount > 0)
            {
                tB->value = (aSum/aCount);
                tB->hasCodValue  = false;
                tB->numCodValues = aCount;
                vSum3 =vSum2*vSum2/aCount;
                if(aCount ==1)
                {
                    tB->sigValue = sqrt(fabs(vSum1-vSum3)/aCount);
                }
                else
                {
                    tB->sigValue = sqrt(fabs(vSum1-vSum3)/(aCount-1));
                }
            }
            else
            {
                std::cout << "Bond DB contains a set with zero member. Bug in the DB"
                          << std::endl;
                exit(1);
            }
        }
        
        
    }
    
    void CodClassify::getCCP4BondAndAngles()
    {
        // Currently CCP4 suite is the only requited thing 
        char * pClibdMon = std::getenv("CLIBD_MON");
        if (pClibdMon !=NULL)
        {
            std::string enerFName(pClibdMon);
            enerFName.append("ener_lib.cif");
            std::ifstream   enerF(enerFName.c_str());
            int nBond = 0;
            int nAng  = 0;
            if (enerF.is_open())
            {
                bool lBond  = false;
                bool lAngle = false;
                bool tOK    = true;       // stop reading if false
                
                std::string tRecord="";
                
                while(!enerF.eof() && tOK)
                {
                    std::getline(enerF, tRecord);
                    tRecord = TrimSpaces(tRecord);
                    
                    std::string tFC=tRecord.substr(0,1);
                    
                    if (tFC.find('#') ==std::string::npos && 
                        tFC.find('.') ==std::string::npos &&  
                        tRecord.size() >0 )
                    {
                        // std::cout << tRecord << std::endl;
                        // std::cout << "C1 " << std::endl;
                        std::vector<std::string> tBuf, tBuf2;
                       
                        if (tRecord.find('#') !=std::string::npos)
                        {
                           StrTokenize(tRecord, tBuf2);
                           if (tBuf2.size() > 1)
                           {
                               StrTokenize(tBuf2[0], tBuf);
                           }
                           else
                           {
                               std::cout<< "What line is that ? " << std::endl 
                                        << std::endl;
                           }
                        }
                        else
                        {
                            StrTokenize(tRecord, tBuf);
                        }
                        if (tBuf.size() > 0)
                        {
                            if (lBond && tBuf.size()==6)
                            {
                                
                                REAL aOrder = StrToOrder(tBuf[2]);
                                ID  aSym    = tBuf[2].substr(0,4);
                                
                                if (aOrder < 0.0)
                                {
                                    std::cout << "Unknown bond order " << tRecord 
                                              << std::endl;
                                    exit(1);
                                }
                                else
                                {
                                    // std::cout << tRecord << std::endl;
                                    if (tBuf[4].size() !=1 )
                                    {
                                        if (tBuf[1].find(".") == std::string::npos)
                                        {
                                            // std::cout << tRecord << std::endl;
                                            ccp4Bonds[tBuf[0]][tBuf[1]]["order"] = aOrder;
                                            ccp4Bonds[tBuf[0]][tBuf[1]]["length"]= StrToReal(tBuf[4]);
                                            
                                            ccp4Bonds[tBuf[1]][tBuf[0]]["order"] = aOrder;
                                            ccp4Bonds[tBuf[1]][tBuf[0]]["length"]= StrToReal(tBuf[4]);
                                            
                                            if (tBuf[5].find('.') == std::string::npos)
                                            {
                                                ccp4Bonds[tBuf[0]][tBuf[1]]["sigValue"] = StrToReal(tBuf[5]);
                                                ccp4Bonds[tBuf[1]][tBuf[0]]["sigValue"] = StrToReal(tBuf[5]);
                                            }
                                            else
                                            {
                                                ccp4Bonds[tBuf[0]][tBuf[1]]["sigValue"] = 0.20;
                                                ccp4Bonds[tBuf[1]][tBuf[0]]["sigValue"] = StrToReal(tBuf[5]);
                                            }
                                            
                                            ccp4BondsA[tBuf[0]][tBuf[1]][aSym]["order"] = aOrder;
                                            ccp4BondsA[tBuf[0]][tBuf[1]][aSym]["length"]= StrToReal(tBuf[4]);
                                            ccp4BondsA[tBuf[1]][tBuf[0]][aSym]["order"] = aOrder;
                                            ccp4BondsA[tBuf[1]][tBuf[0]][aSym]["length"]= StrToReal(tBuf[4]);
                                            
                                            
                                            if (tBuf[5].find('.') == std::string::npos)
                                            {
                                                ccp4BondsA[tBuf[0]][tBuf[1]][aSym]["sigValue"] = StrToReal(tBuf[5]);
                                                ccp4BondsA[tBuf[1]][tBuf[0]][aSym]["sigValue"] = StrToReal(tBuf[5]);
                                            }
                                            else
                                            {
                                                ccp4BondsA[tBuf[0]][tBuf[1]][aSym]["sigValue"] = 0.20;
                                                ccp4BondsA[tBuf[1]][tBuf[0]][aSym]["sigValue"] = 0.20;
                                            }
                                            
                                            nBond++;
                                        }
                                        else
                                        {
                                            ccp4Bonds[tBuf[0]][tBuf[1]]["order"] = aOrder;
                                            ccp4Bonds[tBuf[0]][tBuf[1]]["length"]= StrToReal(tBuf[4]);
                                            if (tBuf[5].find('.') == std::string::npos)
                                            {
                                                ccp4Bonds[tBuf[0]][tBuf[1]]["sigValue"] = StrToReal(tBuf[5]);
                                            }
                                            else
                                            {
                                                ccp4Bonds[tBuf[0]][tBuf[1]]["sigValue"] = 0.20; 
                                            }
                                            
                                            ccp4BondsA[tBuf[0]][tBuf[1]][aSym]["order"] = aOrder;
                                            ccp4BondsA[tBuf[0]][tBuf[1]][aSym]["length"]= StrToReal(tBuf[4]);
                                            if (tBuf[5].find('.') == std::string::npos)
                                            {
                                                ccp4BondsA[tBuf[0]][tBuf[1]][aSym]["sigValue"] = StrToReal(tBuf[5]);
                                            }
                                            else
                                            {
                                                ccp4BondsA[tBuf[0]][tBuf[1]][aSym]["sigValue"] = 0.20; 
                                            }
                                            
                                            nBond++;
                                        }
                                    }
                                }
                            }
                            else if (lAngle && tBuf.size()==5)
                            {
                                ccp4Angles[tBuf[1]][tBuf[0]][tBuf[2]] =  StrToReal(tBuf[4]);
                                nAng++;
                            }
                            else if (tRecord.find("_lib_bond.value_esd") !=std::string::npos)
                            {
                                lBond = true;
                            }
                            else if (tRecord.find("loop_") !=std::string::npos && lBond)
                            {
                                lBond = false;
                            }
                            else if (tRecord.find("_lib_angle.value") !=std::string::npos)
                            {
                                lAngle = true;
                            }
                            else if (tRecord.find("loop_") !=std::string::npos && lAngle)
                            {
                                tOK = false;
                                lAngle = false;
                            }
                        }
                        
                    }
                }
                
                enerF.close();
                
            }
            else
            {
                std::cout << enerFName << " can not be open for reading "
                          << std::endl;
                
                exit(1);
            }
            
          
            
            /*
            std::cout << "Number of ccp4 bonds " << nBond << std::endl;
      
            if (ccp4Bonds.size() > 0)
            {
                std::cout << "The following are the bonds between atoms of different CCP4 types"
                          << std::endl;
                for (std::map<ID, std::map<ID, std::map<ID, REAL> > >::iterator iBA1=ccp4Bonds.begin();
                         iBA1 != ccp4Bonds.end(); iBA1++)
                {
                    for (std::map<ID, std::map<ID, REAL> >::iterator iBA2=iBA1->second.begin();
                            iBA2 !=iBA1->second.end(); iBA2++)
                    {
                        std::cout << iBA1->first << " ----- " << iBA2->first  << std::endl
                                  << " order  : " << iBA2->second["order"]    << std::endl
                                  << " length : " << iBA2->second["length"]   << std::endl
                                  << " stdVar : " << iBA2->second["sigValue"] << std::endl;
                    }  
                }   
            }
            
            std::cout << "Number of ccp4 Angles " << nAng << std::endl;
            if (ccp4Angles.size() > 0)
            {
                std::cout << "The following are the angles between atoms of different CCP4 types "
                          << std::endl;
            
                for (std::map<ID, std::map<ID, std::map<ID, REAL> > >::iterator iAn=ccp4Angles.begin();
                    iAn !=ccp4Angles.end(); iAn++)
                {
                    for (std::map<ID, std::map<ID, REAL> >::iterator iAn1=iAn->second.begin();
                             iAn1 !=iAn->second.end(); iAn1++)
                    {
                        for (std::map<ID, REAL>::iterator iAn2=iAn1->second.begin(); 
                                 iAn2 !=iAn1->second.end(); iAn2++)
                        {
                            std::cout << "Center atom " << iAn1->first << " and atoms "
                                      << iAn->first << " and " << iAn2->first << std::endl 
                                      << iAn2->second << std::endl;
                        }
                    }
                }
            }
            */
        }
        else
        {
            std::cerr << "You need to setup CCP4 suite first " << std::endl;
            exit(1);
        }
        
       
        
        
    }
    
    void CodClassify::setupTargetBonds()
    {    
        
        groupCodOrgBonds();
        
        groupCodMetBonds();
        searchCodBonds();
        constrBondSigmas();
    }
    
    void CodClassify::setupTargetBonds2()
    {    
        
        groupCodOrgBonds2();
        
        groupCodMetBonds();
        searchCodBonds();
        constrBondSigmas();
    }
    /*
    void CodClassify::setupTargetBondsUsingSqlite()
    {
        groupCodMetBonds();
        searchCodBondsUsingSqlite3();  
        constrBondSigmas();
    }
    */
    /*
    void CodClassify::searchCodBondsUsingSqlite3()
    {
        
        sqlite3 * combDB;
        
        int rC=0;   
        
        //std::string clibMonDir(std::getenv("CLIBD_MON"));
        //std::string combDBName = clibMonDir + "allOrg.db";
        //std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
        std::string combDBName = clibMonDir + "/tables/allOrg.db";
        
        rC=sqlite3_open(combDBName.c_str(), &combDB);
        if (rC)
        {
            std::cout << "Can't open the database: " << combDBName <<std::endl;
            sqlite3_close(combDB);
            exit(1);
        }
        
        int nFind =0;
        for (std::vector<BondDict>::iterator iGB = allBonds.begin();
                iGB != allBonds.end(); iGB++)
        {
            std::cout << "Bond between "<< iGB->atoms[0] << " and "
                    << iGB->atoms[1] << std::endl;
            
            bool mB = false;
            for (std::map<ID, int>::iterator iKey=iGB->fullAtoms.begin();
                    iKey !=iGB->fullAtoms.end(); iKey++)
            {
                if (allAtoms[iKey->second].isMetal)
                {
                    mB = true;
                }
            }
            
          
            if (mB)
            {
                std::cout << " search a metal-related bond " << std::endl;
                searchCodMetBonds(iGB);
            }
            else
            {
                std::cout << " search for a organic bond " << std::endl;
                searchOneOrgBondFromCodUsingSqlite(combDB, iGB);
            }
            std::cout << "The final target bond value is " 
                      << iGB->value << std::endl;
            if (iGB->value > 0.3)
            {
                nFind++;
            }
        }   
        std::cout << "Number of Bonds to be found " << (int)allBonds.size() << std::endl
                  << "Number of Bonds found " << nFind << std::endl;
        
        
        sqlite3_close(combDB);
         
    }
    */
   
    /*
    void CodClassify::searchOneOrgBondFromCodUsingSqlite(sqlite3 * tCombDB,  
                                                         std::vector<BondDict>::iterator tB)
    {
        
        std::vector<int> tPair;
        for (std::map<ID, int>::iterator iA=tB->fullAtoms.begin();
                iA !=tB->fullAtoms.end(); iA++)
        {
            tPair.push_back(iA->second);
        }
        int ha1, ha2;
        ID a1NB2, a2NB2, a1NB, a2NB, a1C, a2C;
            
        if((int) allAtoms[tPair[0]].hashingValue < 
                   (int) allAtoms[tPair[1]].hashingValue )
        {
            ha1  =allAtoms[tPair[0]].hashingValue;
            ha2  =allAtoms[tPair[1]].hashingValue;
            a1NB2= allAtoms[tPair[0]].codNB2Symb;
            a2NB2= allAtoms[tPair[1]].codNB2Symb;
            a1NB = allAtoms[tPair[0]].codNBSymb;
            a2NB = allAtoms[tPair[1]].codNBSymb;
            a1C  = allAtoms[tPair[0]].codClass;
            a2C  = allAtoms[tPair[1]].codClass;
        }
        else if ((int) allAtoms[tPair[0]].hashingValue == 
                 (int) allAtoms[tPair[1]].hashingValue)
        {
            if((int)allAtoms[tPair[0]].codClass.size() <= 
               (int)allAtoms[tPair[1]].codClass.size())
            {
                ha1  =allAtoms[tPair[0]].hashingValue;
                ha2  =allAtoms[tPair[1]].hashingValue;
                a1NB2= allAtoms[tPair[0]].codNB2Symb;
                a2NB2= allAtoms[tPair[1]].codNB2Symb;
                a1NB = allAtoms[tPair[0]].codNBSymb;
                a2NB = allAtoms[tPair[1]].codNBSymb;
                a1C  = allAtoms[tPair[0]].codClass;
                a2C  = allAtoms[tPair[1]].codClass;
            }
            else
            {
                ha1  =allAtoms[tPair[1]].hashingValue;
                ha2  =allAtoms[tPair[0]].hashingValue;
                a1NB2= allAtoms[tPair[1]].codNB2Symb;
                a2NB2= allAtoms[tPair[0]].codNB2Symb;
                a1NB = allAtoms[tPair[1]].codNBSymb;
                a2NB = allAtoms[tPair[0]].codNBSymb;
                a1C  = allAtoms[tPair[1]].codClass;
                a2C  = allAtoms[tPair[0]].codClass;
            }
        }
        else
        {
            ha1  =allAtoms[tPair[1]].hashingValue;
            ha2  =allAtoms[tPair[0]].hashingValue;
            a1NB2= allAtoms[tPair[1]].codNB2Symb;
            a2NB2= allAtoms[tPair[0]].codNB2Symb;
            a1NB = allAtoms[tPair[1]].codNBSymb;
            a2NB = allAtoms[tPair[0]].codNBSymb;
            a1C  = allAtoms[tPair[1]].codClass;
            a2C  = allAtoms[tPair[0]].codClass;
        }
        
        std::map<ID, ID>  propNB;
        propNB["a1NB2"] = a1NB2;
        propNB["a2NB2"] = a2NB2;
        propNB["a1NB"]  = a1NB;
        propNB["a2NB"]  = a2NB;
        
        std::map<ID, ID> propHash;
        propHash["ha1"] = IntToStr(ha1);
        propHash["ha2"] = IntToStr(ha2);
        
        int dLev = 0;
            
        std::cout << "for target bond of atoms " <<tB->atoms[0] << " and "
                  << tB->atoms[1] <<std::endl;
                       
        for (std::map<ID, int>::iterator iAt=tB->fullAtoms.begin();
               iAt != tB->fullAtoms.end(); iAt++)
        {
            std::cout << "chemType : " << allAtoms[iAt->second].chemType << " COD classes : " 
                      << allAtoms[iAt->second].codClass << std::endl;
        }
                
        std::cout << "ha1 "    << ha1 << " ha2 " << ha2 << std::endl
                  << " a1NB2 " << a1NB2 << " a2NB2 " << a2NB2  << std::endl
                  << " a1NB "  << a1NB  << " a2NB " << a2NB << std::endl
                  << " a1C "   << a1C   << " a2C "  << a2C << std::endl;
        
        // Now query the databases
        
        // 1. search for exact match of atom classes.
        std::string qQue = "SELECT * from bonds WHERE atom1=\'" 
                           + a1C + "\' and atom2=\'" + a2C + "\';"; // one query statement 
        // std::cout << "1st qQue is " << qQue << std::endl;
        
        
        std::vector<std::vector<std::string> > qResults;
        
        sqlite3Query(tCombDB, qQue.c_str(), qResults);
        
        if ((int)qResults.size() !=0)
        {
            // find exact match
            int nVals= StrToInt(qResults[0][10]);
            if (nVals >=4)
            {
                tB->valueST      = StrToReal(qResults[0][8]);
                tB->value        = tB->valueST;
                tB->sigValueST   = StrToReal(qResults[0][9]);
                tB->sigValue     = tB->sigValueST;
                tB->numCodValues = nVals;
            }
            else
            {
               
                dLev = 5;
                searchOneOrgBondUsingSqliteL(tCombDB, tB, dLev, propNB, propHash); 
            }
        }
        else
        {
            dLev = 5;
            searchOneOrgBondUsingSqliteL(tCombDB, tB, dLev, propNB, propHash);
        }  
        
        
    }
    */
    
    // a recursive function for searching 
    /*
    bool CodClassify::searchOneOrgBondUsingSqliteL(sqlite3* tCombDB,  
                                                   std::vector<BondDict>::iterator tB, 
                                                   int  & tLev,
                                                   std::map<ID, ID>  tPropNB,
                                                   std::map<ID, ID> tPropHash)
    {
        
        
        bool iFind = false;
        std::string tQue;
        std::vector<std::vector<std::string> > qResults;
        
        // set different query strings for different levels of searching
        if (tLev >0)
        {
            setQueStr(tQue, tLev, tPropNB, tPropHash);
            sqlite3Query(tCombDB, tQue.c_str(), qResults);
            if ((int)qResults.size() !=0)
            { 
                setOneBondByMean(qResults, tB);
                iFind = true;
            }
            else
            {
               tLev -=1;
               iFind = searchOneOrgBondUsingSqliteL(tCombDB,  tB, tLev,
                                                    tPropNB, tPropHash);
            }
        }
          
        return iFind;
         
        

    }
    */
    void CodClassify::setQueStr(std::string & tQue, int tLev,
                                std::map<ID, ID>  tPropNB,
                                std::map<ID, ID> tPropHash)
    {
        if (tLev==5)
        {
            tQue = "SELECT * from bonds WHERE hash1=" 
                    + tPropHash["ha1"] + " and hash2=" + tPropHash["ha2"]
                    + " and NeighB1N=\'" + tPropNB["a1NB2"] + "\' and NeighB2N=\'" + tPropNB["a2NB2"]
                    + "\' and NeighB1C=\'" + tPropNB["a1NB"] +"\' and NeighB2C=\'" + tPropNB["a2NB"]
                    + "\';";
        }
        else if(tLev==4)
        {
            tQue = "SELECT * from bonds WHERE hash1=" 
                    + tPropHash["ha1"] + " and hash2=" + tPropHash["ha2"]
                    + " and NeighB1N=\'" + tPropNB["a1NB2"] + "\' and NeighB2N=\'" + tPropNB["a2NB2"]
                    + "\' and NeighB1C=\'" + tPropNB["a1NB"] + "\';";
        }
        else if(tLev==3)
        {
            tQue = "SELECT * from bonds WHERE hash1=" 
                    + tPropHash["ha1"] + " and hash2=" + tPropHash["ha2"]
                    + " and NeighB1N=\'" + tPropNB["a1NB2"] + "\' and NeighB2N=\'" + tPropNB["a2NB2"]
                    + "\';";
        }
        else if(tLev==2)
        {
            tQue = "SELECT * from bonds WHERE hash1=" 
                    + tPropHash["ha1"] + " and hash2=" + tPropHash["ha2"]
                    + " and NeighB1N=\'" + tPropNB["a1NB2"] 
                    + "\';";
        }
        else if(tLev==1)
        {
            tQue = "SELECT * from bonds WHERE hash1=" 
                    + tPropHash["ha1"] + " and hash2=" + tPropHash["ha2"]
                    + ";";
        }
        
    }
    
    void CodClassify::setOneBondByMean(std::vector<std::vector<std::string> >& tResults, 
                                 std::vector<BondDict>::iterator tB)
    {
        int nBonds=0;
        REAL sumVals=0.0;
        for (int i=0; i < (int)tResults.size(); i++)
        {
            int n1       = StrToInt(tResults[i][10]); 
            nBonds  += n1;
            sumVals +=(StrToReal(tResults[i][8])*n1);
        }
        
        if (nBonds >0)
        {
            REAL tMean = sumVals/nBonds;
            tB->valueST = tMean;
            tB->value   = tMean;
            REAL tCur  = 0.0;
            for (int i=0; i < (int)tResults.size(); i++)
            {
                int n1  = StrToInt(tResults[i][10]);
                REAL m1 = StrToReal(tResults[i][8]);
                for (int j=0; j < n1; j++)
                {
                    tCur +=((tMean-m1)*(tMean-m1));
                }
            }
            
            if (nBonds ==1)
            {
                tB->sigValue = sqrt(tCur);
            }
            else
            {
                tB->sigValueST = sqrt(tCur/(nBonds-1));
            }
            tB->sigValue = tB->sigValueST;
        }
        
    }
    
    void CodClassify::setOneBondByDist(std::vector<std::vector<std::string> >& tResults, 
                                       std::vector<BondDict>::iterator tB)
    {
    }
    /*
    void CodClassify::sqlite3Query(sqlite3 *      tDB,
                                   SqliteStatment tQue,
                                   std::vector<std::vector<std::string> >  & tResults)
    {
        sqlite3_stmt *statement;
        
        if(sqlite3_prepare_v2(tDB, tQue, -1, &statement, 0) == SQLITE_OK)
        {
            int cols = sqlite3_column_count(statement);
            int result = 0;
            while(true)
            {
                result = sqlite3_step(statement);
                
                if(result == SQLITE_ROW)
                {
                    std::vector<std::string> values;
                    for(int col = 0; col < cols; col++)
                    {
                        std::string  val;
                        char * ptr = (char*)sqlite3_column_text(statement, col);
                        if(ptr)
                        {
                            val = ptr;
                        }
                        else 
                        {
                            val = ""; 
                        }
                        values.push_back(val);
                    }
                    tResults.push_back(values);
                }
                else
                {
                    break;  
                }
            }    
            
            sqlite3_finalize(statement);
        }     
        
        std::string error = sqlite3_errmsg(tDB);
        if(error != "not an error")
        {
            std::cout << tQue << " " << error << std::endl;
        }   
    }
    */
    void CodClassify::constrBondSigmas()
    {
        // a temporary functions for constrain the sigmas of the bonds 
        // in-between 0.02 - 0.01 (This is just REFMAC optimization, 
        // should be cancel when the new geometrical optimizer is done
        for (std::vector<BondDict>::iterator iBo = allBonds.begin();
                iBo != allBonds.end(); iBo++)
        {
            if (iBo->sigValue > 0.02)
            {
                iBo->sigValue = 0.02;
            }
            else if (iBo->sigValue <0.01)
            {
                iBo->sigValue = 0.01;
            }          
        }
    }
    
    void CodClassify::addAtomClassToBonds(std::vector<BondDict> & tBonds)
    {
        for (std::vector<BondDict>::iterator iBo = tBonds.begin();
                iBo != tBonds.end(); iBo++)
        {
            if ((int)iBo->atomsCodClasses.size() ==0)
            {
                iBo->atomsCodClasses.push_back("");
                iBo->atomsCodClasses.push_back("");
            }
            
            for (std::vector<AtomDict>::iterator iAt = allAtoms.begin();
                    iAt != allAtoms.end(); iAt++)
            {
                if (iBo->atoms[0].compare(iAt->id)==0)
                {
                    iBo->atomsCodClasses[0] = iAt->codClass;
                }
                else if (iBo->atoms[1].compare(iAt->id)==0)
                {
                    iBo->atomsCodClasses[1] = iAt->codClass;
                }
            }
        }
        
        // check
        std::cout <<std::endl;
        
        for (std::vector<BondDict>::iterator iBo = tBonds.begin();
                iBo != tBonds.end(); iBo++)
        {
            std::cout<< "Bond : " << iBo->seriNum << std::endl;
            std::cout <<"Its reside : " << iBo->resName << std::endl;
            std::cout << "Its atomic element IDs are " << std::endl;
            std::cout << iBo->atoms[0] << "  and  " << iBo->atoms[1] << std::endl;
            std::cout << "Its COD IDs are " << std::endl
                      << iBo->atomsCodClasses[0] << std::endl
                      << iBo->atomsCodClasses[1] << std::endl;
            std::cout <<std::endl;                     
        }    
    }
    
    
    // Angles related
    
    void CodClassify::initTargetAngles()
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
                  << " in the system " << std::endl;
        
        std::cout << "These angles are : " << std::endl;
        
        for (std::map<int, std::vector<std::vector<int> > >::iterator iAG 
               =allAnglesIdxs.begin(); iAG != allAnglesIdxs.end(); iAG++)
        {
            std::cout << "There are " << (int)iAG->second.size() 
            << " angles centered on atom " << allAtoms[iAG->first].id << std::endl
            << "They are " << std::endl;
                    
            for(std::vector<std::vector<int> >::iterator iAN =iAG->second.begin();
                    iAN != iAG->second.end(); iAN++)
            {
                std::cout << "atoms: " << allAtoms[iAG->first].id <<", COD code  " 
                          << allAtoms[iAG->first].codClass << std::endl;
                for (std::vector<int>::iterator iAt = iAN->begin(); 
                        iAt !=iAN->end(); iAt++)
                {
                   std::cout << allAtoms[*iAt].id << ", COD code " 
                          << allAtoms[*iAt].codClass << std::endl;
                }
                std::cout<<std::endl;
            }
        }
        
    }
    
    void CodClassify::setDefaultOrgAngle()
    {
        // The key are related to sp hybridization,
        // the values are the default angles
        DefaultOrgAngles[1] = 180.000;
        DefaultOrgAngles[2] = 120.000;
        DefaultOrgAngles[3] = 109.471; 
    }
    
    void CodClassify::setDefaultCoordGeos()
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
        //std::string metDefCoordGeoFileName = clibMonDir + "allMetalDefCoordGeos.table";
        //std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
        std::string metDefCoordGeoFileName = libmolTabDir  + "/allMetalDefCoordGeos.table";
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
    
    void CodClassify::groupCodOrgAngles()
    {
        try
        {
            //std::string clibMonDir(std::getenv("CLIBD_MON"));
            //std::string codAngleFileName = clibMonDir + "allOrgAngles.table";
            //std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
            std::string codAngleFileName = libmolTabDir  + "/allOrgAngles.table";
            std::ifstream codAngleFile(codAngleFileName.c_str());
            
            if(codAngleFile.is_open())
            {
                std::string tRecord="";
                int nLine =0;
        
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
                                =nLine;
                        
                        // allDictAnglesIdx2[ha1][ha2][ha3][tBuf[9]][tBuf[10]][tBuf[11]].push_back(nLine);

                        AngleDict  aAngle;
                        aAngle.seriNum    = nLine +1;
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
                        
                        nLine +=1;
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
    
    void CodClassify::setOrgAngleHeadHashList()
    {
        // std::string clibMonDir(std::getenv("CLIBD_MON"));
        // std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
        for (std::vector<AngleDict>::iterator iAn = allAngles.begin();
                iAn != allAngles.end(); iAn++)
        {
            std::string fRoot = libmolTabDir  + + "/allOrgAngleTables/";
            int iSmall = 2*wSize;
            for (std::vector<int>::iterator iA=iAn->atoms.begin();
                    iA != iAn->atoms.end(); iA++)
            {
                if (allAtoms[*iA].hashingValue < iSmall)
                {
                    iSmall = allAtoms[*iA].hashingValue;
                }
            }
            std::map<int, ID>::iterator iFind = codOrgAngleFiles.find(iSmall);
            if (iFind ==codOrgAngleFiles.end())
            {
                codOrgAngleFiles[iSmall] = fRoot + IntToStr(iSmall) + ".table";
            }
        }
        /*
        std::cout << "Angles are in the following files :" << std::endl;
        for (std::map<int, ID>::iterator iAF = codOrgAngleFiles.begin();
                iAF !=codOrgAngleFiles.end(); iAF++)
        {
            std::cout << iAF->second << std::endl;
        } 
         */        
    }
 
    
    void CodClassify::setOrgAngleHeadHashList2()
    {
        std::map<ID, std::map<ID, std::map<ID, ID> > > allAngIdx;
        
        //std::string clibMonDir(std::getenv("CLIBD_MON"));
        //std::string fRoot = clibMonDir + "allOrgAngleTables/";
        // std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
        std::string fRoot = libmolTabDir  + "/allOrgAngleTables/";
        std::string fIdx  = fRoot + "angle_idx.table";
        std::ifstream codAngleIdxFile(fIdx.c_str());
        if (codAngleIdxFile.is_open())
        {
            std::string tRecord="";
            while(!codAngleIdxFile.eof())
            {
                std::getline(codAngleIdxFile, tRecord);
                if (tRecord.size() !=0)
                {
                    tRecord = TrimSpaces(tRecord);
                    std::vector<std::string> tBuf;
                    StrTokenize(tRecord, tBuf);
                    if ((int)tBuf.size() ==4)
                    {
                        allAngIdx[tBuf[0]][tBuf[1]][tBuf[2]] = tBuf[3];
                    }
                }
            }
            codAngleIdxFile.close();
        }
        
        
        for (std::vector<AngleDict>::iterator iAn = allAngles.begin();
                iAn != allAngles.end(); iAn++)
        {
            //std::string fRoot = "/Users/flong/COD/New_EXP/Current/derivedData/allOrgBondTables
            
            ID ha0, ha1, ha2;
            ha0 = IntToStr(allAtoms[iAn->atoms[0]].hashingValue);
            if (allAtoms[iAn->atoms[1]].hashingValue <= allAtoms[iAn->atoms[2]].hashingValue)
            {
                ha1 = IntToStr(allAtoms[iAn->atoms[1]].hashingValue);
                ha2 = IntToStr(allAtoms[iAn->atoms[2]].hashingValue);
            }
            else
            {
                ha1 = IntToStr(allAtoms[iAn->atoms[2]].hashingValue);
                ha2 = IntToStr(allAtoms[iAn->atoms[1]].hashingValue);
            }
                
            ID haNum;
            if (allAngIdx.find(ha0) != allAngIdx.end())
            {
                if(allAngIdx[ha0].find(ha1) !=allAngIdx[ha0].end())
                {
                    if (allAngIdx[ha0][ha1].find(ha2) !=allAngIdx[ha0][ha1].end())
                    {
                        haNum = allAngIdx[ha0][ha1][ha2];
                    
                        std::map<ID, ID>::iterator iFind = codOrgAngleFiles2.find(haNum);
                        // std::cout << "ha0 " << ha0 << " ha1 " << ha1 << " ha2 " << ha2 << " haNum " << haNum 
                        //          << std::endl;
                        if (iFind ==codOrgAngleFiles2.end())
                        {
                            codOrgAngleFiles2[haNum] = fRoot + haNum + ".table";
                        }
                    }
                }
            }
        }
        /*
        std::cout << "Angles are in the following files :" << std::endl;
        for (std::map<ID, ID>::iterator iAF = codOrgAngleFiles2.begin();
                iAF !=codOrgAngleFiles2.end(); iAF++)
        {
            std::cout << iAF->second << std::endl;
        } 
         */
        
    }
 
    void CodClassify::setOrgAngleHeadHashList22()
    {
        std::map<ID, std::map<ID, std::map<ID, ID> > > allAngIdx;
        
        //std::string clibMonDir(std::getenv("CLIBD_MON"));
        //std::string fRoot = clibMonDir + "allOrgAngleTables/";
        // std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
        std::string fRoot = libmolTabDir  + "/allOrgAngleTables/";
        // std::string fRoot = "/Applications/ccp4-6.5/share/acedrg/tables_DB3/allOrgAngleTables/";
        std::string fIdx  = fRoot + "angle_idx.table";
        std::ifstream codAngleIdxFile(fIdx.c_str());
        if (codAngleIdxFile.is_open())
        {
            std::string tRecord="";
            while(!codAngleIdxFile.eof())
            {
                std::getline(codAngleIdxFile, tRecord);
                if (tRecord.size() !=0)
                {
                    tRecord = TrimSpaces(tRecord);
                    std::vector<std::string> tBuf;
                    StrTokenize(tRecord, tBuf);
                    if ((int)tBuf.size() ==4)
                    {
                        allAngIdx[tBuf[0]][tBuf[1]][tBuf[2]] = tBuf[3];
                    }
                }
            }
            codAngleIdxFile.close();
        }
        
        
        for (std::vector<AngleDict>::iterator iAn = allAngles.begin();
                iAn != allAngles.end(); iAn++)
        {
            //std::string fRoot = "/Users/flong/COD/New_EXP/Current/derivedData/allOrgBondTables
            
            ID ha0, ha1, ha2;
            ha0 = IntToStr(allAtoms[iAn->atoms[0]].hashingValue);
            if (allAtoms[iAn->atoms[1]].hashingValue <= allAtoms[iAn->atoms[2]].hashingValue)
            {
                ha1 = IntToStr(allAtoms[iAn->atoms[1]].hashingValue);
                ha2 = IntToStr(allAtoms[iAn->atoms[2]].hashingValue);
            }
            else
            {
                ha1 = IntToStr(allAtoms[iAn->atoms[2]].hashingValue);
                ha2 = IntToStr(allAtoms[iAn->atoms[1]].hashingValue);
            }
                
            ID haNum;
            if (allAngIdx.find(ha0) != allAngIdx.end())
            {
                if(allAngIdx[ha0].find(ha1) !=allAngIdx[ha0].end())
                {
                    if (allAngIdx[ha0][ha1].find(ha2) !=allAngIdx[ha0][ha1].end())
                    {
                        haNum = allAngIdx[ha0][ha1][ha2];
                    
                        std::map<ID, ID>::iterator iFind = codOrgAngleFiles2.find(haNum);
                        // std::cout << "ha0 " << ha0 << " ha1 " << ha1 << " ha2 " << ha2 << " haNum " << haNum 
                        //          << std::endl;
                        if (iFind ==codOrgAngleFiles2.end())
                        {
                            codOrgAngleFiles2[haNum] = fRoot + haNum + ".table";
                        }
                    }
                }
            }
        }
        
        std::cout << "Angles are in the following files :" << std::endl;
        for (std::map<ID, ID>::iterator iAF = codOrgAngleFiles2.begin();
                iAF !=codOrgAngleFiles2.end(); iAF++)
        {
            std::cout << iAF->second << std::endl;
        } 
        
        
    }
 
    
    void CodClassify::groupCodOrgAngles2()
    {
        setOrgAngleHeadHashList2();
        time_t tStart, tEnd;
        std::time (&tStart);
        std::cout << "Clustering all COD angles started at " << std::ctime(&tStart);
        int nLine =0;
        for (std::map<ID, ID>::iterator iAF=codOrgAngleFiles2.begin();
                    iAF !=codOrgAngleFiles2.end(); iAF++)
        {   
            std::ifstream codAngleFile(iAF->second.c_str());
            if(codAngleFile.is_open())
            {
                // should be something like std::string tNewCodBondFileName(clibMonDir + "/list/bonds.txt");
                // std::string clibMonDir(std::getenv("CLIBD_MON"));
                // std::string codAngleFileName = clibMonDir + "/angles_remain.table";
                // std::ifstream codAngleFile(codAngleFileName.c_str());
                // std::ifstream codAngleFile("/Users/flong/COD/New_EXP/Current/derivedData/allOrgAngles.table");
            
                std::string tRecord="";
                // int nLine =0;
                
                while(!codAngleFile.eof())
                {
                    std::getline(codAngleFile, tRecord);
                    
                    tRecord = TrimSpaces(tRecord);
                    std::vector<std::string> tBuf;
                    StrTokenize(tRecord, tBuf);
                    
                    if ((int)tBuf.size() ==24)
                    {
                        int ha1, ha2, ha3;
                        
                        ha1 = StrToInt(tBuf[0]);
                        ha2 = StrToInt(tBuf[1]);
                        ha3 = StrToInt(tBuf[2]);
                        
                        //std::cout << ha1 << std::endl;
                        /*
                        for (int i=3; i < (int)tBuf.size(); i++)
                        {
                            tBuf[i] = TrimSpaces(tBuf[i]);
                        }
                        */
                        
                        allDictAnglesIdx[ha1][ha2][ha3][tBuf[3]][tBuf[4]][tBuf[5]][tBuf[6]][tBuf[7]][tBuf[8]][tBuf[9]][tBuf[10]][tBuf[11]]
                                =nLine;
                        
                        aValueSet aAngS3;
                        aAngS3.value        = StrToReal(tBuf[21]);
                        aAngS3.sigValue     = StrToReal(tBuf[22]);
                        aAngS3.numCodValues = StrToInt(tBuf[23]);
                        if (allDictAnglesIdx3[ha1][ha2][ha3].size()==0)
                        {
                            allDictAnglesIdx3[ha1][ha2][ha3].push_back(aAngS3);
                        }
                        
                        aValueSet aAngS2;
                        aAngS2.value        = StrToReal(tBuf[18]);
                        aAngS2.sigValue     = StrToReal(tBuf[19]);
                        aAngS2.numCodValues = StrToInt(tBuf[20]);
                        if (allDictAnglesIdx2[ha1][ha2][ha3][tBuf[3]][tBuf[4]][tBuf[5]].size()==0)
                        {
                            allDictAnglesIdx2[ha1][ha2][ha3][tBuf[3]][tBuf[4]][tBuf[5]].push_back(aAngS2);
                        }
                        
                        aValueSet aAngS1;
                        aAngS1.value        = StrToReal(tBuf[15]);
                        aAngS1.sigValue     = StrToReal(tBuf[16]);
                        aAngS1.numCodValues = StrToInt(tBuf[17]);
                        if (allDictAnglesIdx1[ha1][ha2][ha3][tBuf[3]][tBuf[4]][tBuf[5]][tBuf[6]][tBuf[7]][tBuf[8]].size()==0)
                        {
                            allDictAnglesIdx1[ha1][ha2][ha3][tBuf[3]][tBuf[4]][tBuf[5]][tBuf[6]][tBuf[7]][tBuf[8]].push_back(aAngS1);
                        }
                        
                        AngleDict  aAngle;
                        aAngle.seriNum    = nLine;
                        aAngle.atomsNB2Rep.push_back(tBuf[3]);
                        aAngle.atomsNB2Rep.push_back(tBuf[4]);
                        aAngle.atomsNB2Rep.push_back(tBuf[5]);
                        aAngle.atomsNBRep.push_back(tBuf[6]);
                        aAngle.atomsNBRep.push_back(tBuf[7]);
                        aAngle.atomsNBRep.push_back(tBuf[8]);
                        
                        aAngle.atomsCodClasses.push_back(tBuf[9]);
                        aAngle.atomsCodClasses.push_back(tBuf[10]);
                        aAngle.atomsCodClasses.push_back(tBuf[11]);
                        
                        aAngle.value        = StrToReal(tBuf[12]);
                        aAngle.sigValue     = StrToReal(tBuf[13]);
                        //aAngle.valueP       = StrToReal(tBuf[15]);
                        //aAngle.sigValueP    = StrToReal(tBuf[16]);
                        if(aAngle.sigValue <0.0001 || aAngle.sigValue > 3.0)
                        {
                            aAngle.sigValue = 3.0;
                        }
                        aAngle.numCodValues  = StrToInt(tBuf[14]);
                        // aAngle.numCodValuesP = StrToInt(tBuf[17]);
 
                        allDictAngles.push_back(aAngle);
                        
                        nLine +=1;
                       
                    }
                     
                }
                
                codAngleFile.close();
            }
            else
            {
                std::cout << iAF->second 
                //<< " can not be open for reading. Check if $CLIBD_MON/allOrgAngleTables contains that file! "
                << " can not be open for reading. Check if "
                << libmolTabDir << "/allOrgAngleTables contains that file! "
                << std::endl;
                exit(1);
            }
        }
        
        std::time(&tEnd);
        std::cout << "Clustering COD org angles finished at " << std::ctime(&tEnd);
        REAL tDiff;
        tDiff = std::difftime(tEnd,tStart);
        std::cout  << "it takes " << std::setprecision(3) <<tDiff 
                   << " seconds to finish group COD angles " << std::endl;
    }
    
    void CodClassify::groupCodOrgAngles22()
    {
        setOrgAngleHeadHashList22();
        time_t tStart, tEnd;
        std::time (&tStart);
        std::cout << "Clustering all COD angles started at " << std::ctime(&tStart);
        int nLine =0;
        for (std::map<ID, ID>::iterator iAF=codOrgAngleFiles2.begin();
                    iAF !=codOrgAngleFiles2.end(); iAF++)
        {   
            std::ifstream codAngleFile(iAF->second.c_str());
            if(codAngleFile.is_open())
            {
                // should be something like std::string tNewCodBondFileName(clibMonDir + "/list/bonds.txt");
                // std::string clibMonDir(std::getenv("CLIBD_MON"));
                // std::string codAngleFileName = clibMonDir + "/angles_remain.table";
                // std::ifstream codAngleFile(codAngleFileName.c_str());
                // std::ifstream codAngleFile("/Users/flong/COD/New_EXP/Current/derivedData/allOrgAngles.table");
            
                std::string tRecord="";
                // int nLine =0;
                
                while(!codAngleFile.eof())
                {
                    std::getline(codAngleFile, tRecord);
                    
                    tRecord = TrimSpaces(tRecord);
                    std::vector<std::string> tBuf;
                    StrTokenize(tRecord, tBuf);
                    
                    if ((int)tBuf.size() ==27)
                    {
                        int ha1, ha2, ha3;
                        
                        ha1 = StrToInt(tBuf[0]);
                        ha2 = StrToInt(tBuf[1]);
                        ha3 = StrToInt(tBuf[2]);
                        
                        //std::cout << ha1 << std::endl;
                        /*
                        for (int i=3; i < (int)tBuf.size(); i++)
                        {
                            tBuf[i] = TrimSpaces(tBuf[i]);
                        }
                        */
                        // exact match  
                        allDictAnglesIdxD[ha1][ha2][ha3][tBuf[3]][tBuf[4]][tBuf[5]][tBuf[6]][tBuf[7]][tBuf[8]]
                                        [tBuf[9]][tBuf[10]][tBuf[11]][tBuf[12]][tBuf[13]][tBuf[14]]
                                       =nLine;
                        nLine +=1;
                        
                        aValueSet aAng;
                        aAng.value        = StrToReal(tBuf[15]);
                        aAng.sigValue     = StrToReal(tBuf[16]);
                        if(aAng.sigValue <0.0001 || aAng.sigValue > 3.0)
                        {
                            aAng.sigValue = 3.0;
                        }
                        aAng.numCodValues  = StrToInt(tBuf[17]);
                        allDictAnglesD.push_back(aAng);
                        
                        
                        // AxC not found 
                        aValueSet aAngS1;
                        aAngS1.value        = StrToReal(tBuf[18]);
                        aAngS1.sigValue     = StrToReal(tBuf[19]);
                        if(aAngS1.sigValue <0.0001 || aAngS1.sigValue > 3.0)
                        {
                            aAngS1.sigValue = 3.0;
                        }
                        aAngS1.numCodValues = StrToInt(tBuf[20]);
                       
                        if (allDictAnglesIdx1D[ha1][ha2][ha3][tBuf[3]][tBuf[4]][tBuf[5]]
                             [tBuf[6]][tBuf[7]][tBuf[8]][tBuf[9]][tBuf[10]][tBuf[11]].size()==0)
                        {
                            allDictAnglesIdx1D[ha1][ha2][ha3][tBuf[3]][tBuf[4]][tBuf[5]]
                              [tBuf[6]][tBuf[7]][tBuf[8]][tBuf[9]][tBuf[10]][tBuf[11]].push_back(aAngS1);
                        } 
                        
                        // AxM not found 
                        aValueSet aAngS2;
                        aAngS2.value        = StrToReal(tBuf[21]);
                        aAngS2.sigValue     = StrToReal(tBuf[22]);
                        if(aAngS2.sigValue <0.0001 || aAngS2.sigValue > 3.0)
                        {
                            aAngS2.sigValue = 3.0;
                        }
                        aAngS2.numCodValues = StrToInt(tBuf[23]);
                       
                        if (allDictAnglesIdx2D[ha1][ha2][ha3][tBuf[3]][tBuf[4]][tBuf[5]]
                             [tBuf[6]][tBuf[7]][tBuf[8]].size()==0)
                        {
                            allDictAnglesIdx2D[ha1][ha2][ha3][tBuf[3]][tBuf[4]][tBuf[5]]
                              [tBuf[6]][tBuf[7]][tBuf[8]].push_back(aAngS2);
                        } 
                        
                        // axNB not found 
                        aValueSet aAngS3;
                        aAngS3.value        = StrToReal(tBuf[24]);
                        aAngS3.sigValue     = StrToReal(tBuf[25]);
                        if(aAngS3.sigValue <0.0001 || aAngS3.sigValue > 3.0)
                        {
                            aAngS3.sigValue = 3.0;
                        }
                        aAngS3.numCodValues = StrToInt(tBuf[26]);
                        if (allDictAnglesIdx3D[ha1][ha2][ha3][tBuf[3]][tBuf[4]][tBuf[5]].size()==0)
                        {
                            allDictAnglesIdx3D[ha1][ha2][ha3][tBuf[3]][tBuf[4]][tBuf[5]].push_back(aAngS3);
                        }
                        
                    }
                     
                }
                
                codAngleFile.close();
            }
            else
            {
                std::cout << iAF->second 
                //<< " can not be open for reading. Check if $CLIBD_MON/allOrgAngleTables contains that file! "
                << " can not be open for reading. Check if "
                << libmolTabDir << "/allOrgAngleTables contains that file! "
                << std::endl;
                exit(1);
            }
        }
        
        std::time(&tEnd);
        std::cout << "Clustering COD org angles finished at " << std::ctime(&tEnd);
        REAL tDiff;
        tDiff = std::difftime(tEnd,tStart);
        std::cout  << "it takes " << std::setprecision(3) <<tDiff 
                   << " seconds to finish group COD angles " << std::endl;
    }
  
    void CodClassify::groupCodMetAngles()
    {
        try
        {
            //std::string clibMonDir(std::getenv("CLIBD_MON"));
            //std::string tAName = clibMonDir + "allMetalCoordGeoAngles.table";
            //std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
            std::string tAName = libmolTabDir  + "/allMetalCoordGeoAngles.table";
            
            std::ifstream codMetAngleFile(tAName.c_str());
            if(codMetAngleFile.is_open())
            {
                std::string tRecord="";
              
                while(!codMetAngleFile.eof())
                {
                    std::getline(codMetAngleFile, tRecord);
                    tRecord = TrimSpaces(tRecord);
                    std::vector<std::string> tBuf;
                    StrTokenize(tRecord, tBuf, ':');
                    
                    if ((int)tBuf.size() == 2)
                    {
                        std::vector<std::string> tBuf1, tBuf2;
                        tBuf[0] = TrimSpaces(tBuf[0]);
                        StrTokenize(tBuf[0], tBuf1, '_');
                        tBuf[1] = TrimSpaces(tBuf[1]);
                        StrTokenize(tBuf[1], tBuf2);
                        
                        for(int i=0; i < (int)tBuf2.size(); i++)
                        {
                            allDictCoordGeoAngs[StrToInt(tBuf1[1])][tBuf1[2]].push_back(StrToReal(tBuf2[i]));
                        }
                    }
                }
                codMetAngleFile.close();
            }
        }
        catch (std::exception & e)
        {
            std::cout << e.what() << std::endl;
        }
    }
    
    void CodClassify::groupCodAnglesWithNonCenteredMetal()
    {
        try
        {
            // should be something like std::string tNewCodBondFileName(clibMonDir + "/list/bonds.txt");
            // std::string clibMonDir(std::getenv("CLIBD_MON"));
            // std::string tNCName =  clibMonDir + "allOrgAnglesWithNonCenteredMetalNB.table";
            // std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
            std::string tNCName =  libmolTabDir  + "/allOrgAnglesWithNonCenteredMetalNB.table";
            
            std::ifstream codNonCenMetAngleFile(tNCName.c_str());
            if(codNonCenMetAngleFile.is_open())
            {
                std::string tRecord="";
                int nLine = (int)allDictAngles.size();
                // std::cout << "non-centered metal starting line " << nLine << std::endl;
                while(!codNonCenMetAngleFile.eof())
                {
                    std::getline(codNonCenMetAngleFile, tRecord);
                    tRecord = TrimSpaces(tRecord);
                    std::vector<std::string> tBuf;
                    StrTokenize(tRecord, tBuf);
                    
                    if ((int)tBuf.size() ==9)
                    {
                        for (int i=0; i < (int)tBuf.size(); i++)
                        {
                            tBuf[i] = TrimSpaces(tBuf[i]);
                        }
                        
                        allDictNonCenMetAnglesIdx[tBuf[0]][tBuf[1]][tBuf[2]][tBuf[3]][tBuf[4]][tBuf[5]] =nLine;
                                                
                        AngleDict  aAngle;
                        aAngle.seriNum    = nLine+1;
                        
                        aAngle.value        = StrToReal(tBuf[6]);
                        aAngle.sigValue     = StrToReal(tBuf[7]);
                        if(aAngle.sigValue <0.0001)
                        {
                            aAngle.sigValue = 3.0;
                        }
                        aAngle.numCodValues = StrToInt(tBuf[8]);
                        
                        allDictAngles.push_back(aAngle);
                        
                        nLine +=1;
                       
                    }

                }
            }
            
        }
        catch (std::exception & e)
        {
            std::cout << e.what() << std::endl;
        }
        
    }
    
    
    void CodClassify::searchCodAngles()
    {
        std::map<int, std::vector<AngleDict> > specialAngs;
        
        for (std::vector<AngleDict>::iterator iA=allAngles.begin();
                iA !=allAngles.end(); iA++)
        { 
            std::cout << "Angle between " << allAtoms[iA->atoms[0]].id 
                              << "(center) and " << allAtoms[iA->atoms[1]].id
                              << " and " << allAtoms[iA->atoms[2]].id << std::endl;
            
            if (allAtoms[iA->atoms[0]].isMetal)
            {
               std::cout << "getIdealCNGeoAngles " << std::endl;
               //searchCodMetAngles(iA);
               getIdealCNGeoAngles(iA);
               std::cout << "angle candidates: " << std::endl;
               for (std::vector<REAL>::iterator iCA=iA->codAngleValues.begin();
                            iCA !=iA->codAngleValues.end(); iCA++)
               {
                   std::cout << *iCA << std::endl;
               }
            }
            else if (!allAtoms[iA->atoms[1]].isMetal && !allAtoms[iA->atoms[2]].isMetal) 
            {
                //std::cout << "searchCodOrgAngles " << std::endl;
                // searchCodOrgAngles(iA
                bool lSpeAng = checkSpeAng(iA);  
                if (lSpeAng)
                {
                    specialAngs[iA->atoms[0]].push_back(*iA);
                }
                else
                {
                    searchCodOrgAngles2(iA);
                }
            }
            else 
            {
                std::cout << "searchCodAnglesWithNonCenteredMetal " << std::endl;
                searchCodAnglesWithNonCenteredMetal(iA);
            }
            
            if (!allAtoms[iA->atoms[0]].isMetal)
            {
                std::cout << "Target angle value : " << iA->value << std::endl;
            }
            
        }
        
        // Deal with the special angles when other angles are done 
        if (specialAngs.size())
        {
            setSpecialAngles(specialAngs);
        }
        
        // Final check all constraints 
        checkAngConstraints();
        
        
    }
    
    void CodClassify::searchCodAngles2()
    {
        std::map<int, std::vector<AngleDict> > specialAngs;
        
        for (std::vector<AngleDict>::iterator iA=allAngles.begin();
                iA !=allAngles.end(); iA++)
        { 
            
            
            if (allAtoms[iA->atoms[0]].isMetal)
            {
               std::cout << "getIdealCNGeoAngles " << std::endl;
               //searchCodMetAngles(iA);
               getIdealCNGeoAngles(iA);
               std::cout << "angle candidates: " << std::endl;
               for (std::vector<REAL>::iterator iCA=iA->codAngleValues.begin();
                            iCA !=iA->codAngleValues.end(); iCA++)
               {
                   std::cout << *iCA << std::endl;
               }
            }
            else if (!allAtoms[iA->atoms[1]].isMetal && !allAtoms[iA->atoms[2]].isMetal) 
            {
                //std::cout << "searchCodOrgAngles " << std::endl;
                // searchCodOrgAngles(iA
                bool lSpeAng = checkSpeAng(iA);  
                if (lSpeAng)
                {
                    specialAngs[iA->atoms[0]].push_back(*iA);
                }
                else
                {
                    searchCodOrgAngles22(iA);
                }
            }
            else 
            {
                std::cout << "searchCodAnglesWithNonCenteredMetal " << std::endl;
                searchCodAnglesWithNonCenteredMetal(iA);
            }
            
            if (!allAtoms[iA->atoms[0]].isMetal)
            {
                std::cout << "Target angle value : " << iA->value << std::endl;
            }
            
        }
        
        // Deal with the special angles when other angles are done 
        if (specialAngs.size())
        {
            setSpecialAngles(specialAngs);
        }
        
        // Final check all constraints 
        checkAngConstraints();
        
        
    }
    
    
    /*
    void CodClassify::searchCodAnglesUsingSqlite()
    {
        
        sqlite3 * combDB;
        
        int rC=0;   
        // Temp name 
        // FileName  combDBName = "/Users/flong/COD/New_EXP/Current/derivedData/SQLite_related/allOrg.db";
        //std::string clibMonDir(std::getenv("CLIBD_MON"));
        //std::string combDBName = clibMonDir + "allOrg.db";
        //std::string clibMonDir(std::getenv("LIBMOL_ROOT"));
        std::string combDBName = clibMonDir + "/tables/allOrg.db";
        
        rC=sqlite3_open(combDBName.c_str(), &combDB);
        if (rC)
        {
            std::cout << "Can't open the database: " << combDBName <<std::endl;
            sqlite3_close(combDB);
            exit(1);
        }

        
        for (std::vector<AngleDict>::iterator iA=allAngles.begin();
                iA !=allAngles.end(); iA++)
        { 
            std::cout << "Angle between " << allAtoms[iA->atoms[0]].id 
                              << "(center) and " << allAtoms[iA->atoms[1]].id
                              << " and " << allAtoms[iA->atoms[2]].id << std::endl;
            
            if (allAtoms[iA->atoms[0]].isMetal)
            {
               std::cout << "getIdealCNGeoAngles " << std::endl;
               //searchCodMetAngles(iA);
               getIdealCNGeoAngles(iA);
               std::cout << "angle candidates: " << std::endl;
               for (std::vector<REAL>::iterator iCA=iA->codAngleValues.begin();
                            iCA !=iA->codAngleValues.end(); iCA++)
               {
                   std::cout << *iCA << std::endl;
               }
            }
            else if (!allAtoms[iA->atoms[1]].isMetal && !allAtoms[iA->atoms[2]].isMetal) 
            {
                //std::cout << "searchCodOrgAngles using sqlite3" << std::endl;
                // searchCodOrgAngles(iA);
                searchOneOrgAngleFromCodUsingSqlite(combDB, iA);
            }
            else 
            {
                //std::cout << "searchCodAnglesWithNonCenteredMetal " << std::endl;
                searchCodAnglesWithNonCenteredMetal(iA);
            }
            
            if (!allAtoms[iA->atoms[0]].isMetal)
            {
                std::cout << "Target angle value : " << iA->value << std::endl;
            }   
        }
        
    }
     * 
     */
    
    void CodClassify::getIdealCNGeoAngles(std::vector<AngleDict>::iterator iAN)
    {
        if (allAtoms[iAN->atoms[0]].isMetal)
        {
            int cn  = allAtoms[iAN->atoms[0]].connAtoms.size();
            ID  ct  = allAtoms[iAN->atoms[0]].chemType;
            ID  geo = allAtoms[iAN->atoms[0]].metalGeo;
            std::cout << "Atom is a " << ct << std::endl;
            std::cout << "Coordination number is " << cn << std::endl;
            std::cout << "coordination geometry is " << geo << std::endl;
            
            if (!geo.empty())
            {
                // the geometry of the metal atom has been set 
                // (either by default or by user input)
                // For those need on-fly calculations  
                if (geo.compare("PENTAGONAL-ANTIPROSMATIC")==0)
                {
                    // std::cout << "run getCN10PentAntiPros" << std::endl;
                    getCN10PentAntiPros(iAN);
                    //std::cout << "angle candidates: " << std::endl;
                    //for (std::vector<REAL>::iterator iCA=iAN->codAngleValues.begin();
                    //        iCA !=iAN->codAngleValues.end(); iCA++)
                    //{
                    //    std::cout << *iCA << std::endl;
                    //}
                }
                else
                {
                    std::cout << "run searchCodMetAngles " << std::endl;
                    searchCodMetAngles(iAN);
                }
                
            }
        }
        else
        {
            std::cout << "The center atom of the angle not a metal atom." << std::endl
                      << " Bug! Check the list of metal element." << std::endl;
            exit(1);
        }
    }
    
    void CodClassify::searchCodMetAngles(std::vector<AngleDict>::iterator iAN)
    {
        int cn = (int)allAtoms[iAN->atoms[0]].connAtoms.size();
        ID  ct = allAtoms[iAN->atoms[0]].chemType;
        ID  geo;
        if (allAtoms[iAN->atoms[0]].metalGeo.empty())
        {
            std::map<ID, std::map<int,ID> >::iterator iFind1 
                  =DefaultCoordGeos2.find(ct);
            if (iFind1 !=DefaultCoordGeos2.end())
            {
                std::map<int,ID>::iterator iFind2
                = DefaultCoordGeos2[ct].find(cn);
                if (iFind2 !=DefaultCoordGeos2[ct].end())
                {
                    geo=DefaultCoordGeos2[ct][cn];
                }
                else
                {
                    geo= DefaultCoordGeos[cn];
                }
            }
            else
            {
                geo = DefaultCoordGeos[cn];
            }
        }
        else
        {
            geo = allAtoms[iAN->atoms[0]].metalGeo;
        }
        
        std::cout << "Element is " << ct << std::endl;
        std::cout << "Coordination number is " << cn << std::endl;
        std::cout << "Geometry is " << geo << std::endl;
        
        std::map<int, std::map<ID, std::vector<REAL> > >::iterator iFind=allDictCoordGeoAngs.find(cn);
        if (iFind != allDictCoordGeoAngs.end())
        {
            for (std::vector<REAL>::iterator iDA=allDictCoordGeoAngs[cn][geo].begin();
                    iDA != allDictCoordGeoAngs[cn][geo].end(); iDA++)
            {
                iAN->codAngleValues.push_back(*iDA);
            }
        }
    }
    
    void CodClassify::getCN10PentAntiPros(std::vector<AngleDict>::iterator iAN)
    {
        // int cn  = allAtoms[iAN->atoms[0]].connAtoms.size();
        // ID  ct  = allAtoms[iAN->atoms[0]].chemType;
        // 1. Get the bond length between atoms in the 5-membered ring
        //    and the bond length between the metal atom and one of atoms
        //    in the five-membered ring
        
        int iNB1=-1, iNB2=-1; 
        for (int i=0; i <(int)allAtoms[iAN->atoms[0]].connAtoms.size();
                i++)
        {
            iNB1 = allAtoms[iAN->atoms[0]].connAtoms[i];
            for (int j=0; j <(int)allAtoms[iNB1].connAtoms.size(); j++)
            {
                // iNB2 should be connected to both metal atom and iNB1
                int tNB2 = allAtoms[iNB1].connAtoms[j];
                std::vector<int>::iterator iFind;
                iFind=std::find(allAtoms[iAN->atoms[0]].connAtoms.begin(),
                                allAtoms[iAN->atoms[0]].connAtoms.end(), tNB2);
                if(iFind !=allAtoms[iAN->atoms[0]].connAtoms.end())
                {
                    iNB2 = tNB2;
                    break;
                }
            }
            if (iNB1 !=-1 && iNB2 != -1)
            {
                break;
            }
        }
        //std::cout << "iNB1 " << iNB1 << std::endl;
        //std::cout << "iNB2 " << iNB2 << std::endl;
        if (iNB1 !=-1 && iNB2 != -1)
        {
            REAL a, b, c, d;
            ID id1, id2, id3;
            id1 = allAtoms[iAN->atoms[0]].id;
            id2 = allAtoms[iNB1].id;
            id3 = allAtoms[iNB2].id;
            a=getOnebond(id2, id3);
            //std::cout << "Bond between " << id2 << " and " << id3 
            //        << " is " << a << std::endl;
            
            b=getOnebond(id1, id2);
            //std::cout << "Bond between " << id1 << " and " << id2 
            //        << " is " << b << std::endl;
            
            if (a >0.0 && b >0.0)
            {
                
                REAL r1 = degreeToRadians(72.0);
                REAL r2 = degreeToRadians(144.0);
                REAL r3 = degreeToRadians(36.0);
                REAL r4 = degreeToRadians(108.0);
                
                REAL x1, y1, x2, y2, x3, y3, x4, y4, z;
                
                c  = a/(2.0*sin(r3));
                d  = b*b-c*c;
                if (d >0)
                {
                    z   = sqrt(d);
                    x1  = c*cos(r1);
                    y1  = c*sin(r1);

                    x2  = c*cos(r2);
                    y2  = c*sin(r2);
  
                    x3  = c*cos(r3);
                    y3  = c*sin(r3);
 
                    x4  = c*cos(r4);
                    y4  = c*sin(r4);
                    
                    // get one set of atomic coordinates under this geometry
                    std::vector<std::vector<REAL> > tAtoms;
                    std::vector<REAL> atom0;
                    atom0.push_back(0.0);
                    atom0.push_back(0.0);
                    atom0.push_back(0.0);
                    tAtoms.push_back(atom0);
                    
                    std::vector<REAL> atom1;
                    atom1.push_back(x1);
                    atom1.push_back(y1);
                    atom1.push_back(z);
                    tAtoms.push_back(atom1);
                    
                    std::vector<REAL> atom2;
                    atom2.push_back(x2);
                    atom2.push_back(y2);
                    atom2.push_back(z);
                    tAtoms.push_back(atom2);
                    
                    std::vector<REAL> atom3;
                    atom3.push_back(x2);
                    atom3.push_back(-y2);
                    atom3.push_back(z);
                    tAtoms.push_back(atom3);
                    
                    std::vector<REAL> atom4;
                    atom4.push_back(x1);
                    atom4.push_back(-y1);
                    atom4.push_back(z);
                    tAtoms.push_back(atom4);
                    
                    std::vector<REAL> atom5;
                    atom5.push_back(c);
                    atom5.push_back(0.0);
                    atom5.push_back(z);
                    tAtoms.push_back(atom5);
          
                    std::vector<REAL> atom6;
                    atom6.push_back(x3);
                    atom6.push_back(y3);
                    atom6.push_back(-z);
                    tAtoms.push_back(atom6);
                   
                    std::vector<REAL> atom7;
                    atom7.push_back(x4);
                    atom7.push_back(y4);
                    atom7.push_back(-z);
                    tAtoms.push_back(atom7);
          
                    std::vector<REAL> atom8;
                    atom8.push_back(-c);
                    atom8.push_back(0.0);
                    atom8.push_back(-z);
                    tAtoms.push_back(atom8);
                   
                    std::vector<REAL> atom9;
                    atom9.push_back(x3);
                    atom9.push_back(-y3);
                    atom9.push_back(-z);
                    tAtoms.push_back(atom9);
                    
                    std::vector<REAL> atom10;
                    atom10.push_back(x4);
                    atom10.push_back(-y4);
                    atom10.push_back(-z);
                    tAtoms.push_back(atom10);                
                    
                    // Get possible angles 
                    std::vector<REAL> existAngs;
                    
                    for (int i=1; i < (int)tAtoms.size(); i++)
                    {
                        
                        std::vector<REAL> v1;
                        for (int k=0; k<3; k++)
                        {
                           v1.push_back(tAtoms[i][k]-tAtoms[0][k]); 
                        }
                        
                        for (int j=i+1; j < (int)tAtoms.size(); j++)
                        {
                            std::vector<REAL> v2;
                            for (int k=0; k<3; k++)
                            {
                                v2.push_back(tAtoms[j][k]-tAtoms[0][k]);
                            }
                            REAL tAng = RadiansToDegree(getAngle2V(v1,v2));
                            bool lExist=false;
                            for (std::vector<REAL>::iterator iA=existAngs.begin();
                                    iA!=existAngs.end(); iA++)
                            {
                                if (fabs(fabs(tAng)-fabs(*iA)) < 2.0)
                                {
                                    lExist=true;
                                    break;
                                }
                            }
                            
                            if(!lExist)
                            {
                                existAngs.push_back(tAng);
                                iAN->codAngleValues.push_back(tAng);
                            }
                        }
                    }
                    
                }
            }
        }
        
    }
    
    REAL CodClassify::getOnebond(ID aID1, ID aID2)
    {
        for (int iB=0; iB < (int)allBonds.size(); iB++)
        {
            if ( (allBonds[iB].atoms[0]==aID1 && allBonds[iB].atoms[1]==aID2)
                  || (allBonds[iB].atoms[1]==aID1 && allBonds[iB].atoms[0]==aID2) )
            {
                return allBonds[iB].value;
            }
        }
        
        return -1.0;
    }
    
    void CodClassify::searchCodOrgAngles(std::vector<AngleDict>::iterator iAN)
    {    
            /*
            codClassToAtomAng(allAtoms[iAN->atoms[0]].codClass, allAtoms[iAN->atoms[0]]);
            codClassToAtomAng(allAtoms[iAN->atoms[1]].codClass, allAtoms[iAN->atoms[1]]);
            codClassToAtomAng(allAtoms[iAN->atoms[2]].codClass, allAtoms[iAN->atoms[2]]);
            */
        
            int ha1, ha2, ha3;
            ID a1NB2, a1NB, a1C, a2NB2, a2NB, a2C,a3NB2, a3NB, a3C;
            ID id2, id3;
            ha1   = allAtoms[iAN->atoms[0]].hashingValue;
            a1NB2 = allAtoms[iAN->atoms[0]].codNB2Symb;
            a1NB  = allAtoms[iAN->atoms[0]].codNBSymb;
            a1C   = allAtoms[iAN->atoms[0]].codClass;
            
            std::cout << "Atom1 " <<  allAtoms[iAN->atoms[0]].id 
                      << " Its Cod class " <<  a1C <<  " and Hashing "
                      << ha1 <<std::endl;
            std::cout << " Its codNBSymb " <<  a1NB << " its codNB2Symb " << a1NB2 << std::endl;
            
            if ((int)allAtoms[iAN->atoms[1]].hashingValue
                    <(int)allAtoms[iAN->atoms[2]].hashingValue )
            {
                ha2  =allAtoms[iAN->atoms[1]].hashingValue;
                ha3  =allAtoms[iAN->atoms[2]].hashingValue;
           
                a2NB2= allAtoms[iAN->atoms[1]].codNB2Symb;
                a3NB2= allAtoms[iAN->atoms[2]].codNB2Symb;
          
                a2NB = allAtoms[iAN->atoms[1]].codNBSymb;
                a3NB = allAtoms[iAN->atoms[2]].codNBSymb;
         
                a2C  = allAtoms[iAN->atoms[1]].codClass;
                a3C  = allAtoms[iAN->atoms[2]].codClass;
                
                id2  = allAtoms[iAN->atoms[1]].id;
                id3  = allAtoms[iAN->atoms[2]].id;
            }
            else if ((int)allAtoms[iAN->atoms[1]].hashingValue
                    ==(int)allAtoms[iAN->atoms[2]].hashingValue)
            {
                if ((int)allAtoms[iAN->atoms[1]].codClass.size() <=
                        (int)allAtoms[iAN->atoms[2]].codClass.size())
                {
                    ha2  =allAtoms[iAN->atoms[1]].hashingValue;
                    ha3  =allAtoms[iAN->atoms[2]].hashingValue;
           
                    a2NB2= allAtoms[iAN->atoms[1]].codNB2Symb;
                    a3NB2= allAtoms[iAN->atoms[2]].codNB2Symb;
          
                    a2NB = allAtoms[iAN->atoms[1]].codNBSymb;
                    a3NB = allAtoms[iAN->atoms[2]].codNBSymb;
         
                    a2C  = allAtoms[iAN->atoms[1]].codClass;
                    a3C  = allAtoms[iAN->atoms[2]].codClass;
                
                    id2  = allAtoms[iAN->atoms[1]].id;
                    id3  = allAtoms[iAN->atoms[2]].id;
                }
                else
                {
                    ha2  =allAtoms[iAN->atoms[2]].hashingValue;
                    ha3  =allAtoms[iAN->atoms[1]].hashingValue;
           
                    a2NB2= allAtoms[iAN->atoms[2]].codNB2Symb;
                    a3NB2= allAtoms[iAN->atoms[1]].codNB2Symb;
          
                    a2NB = allAtoms[iAN->atoms[2]].codNBSymb;
                    a3NB = allAtoms[iAN->atoms[1]].codNBSymb;
         
                    a2C  = allAtoms[iAN->atoms[2]].codClass;
                    a3C  = allAtoms[iAN->atoms[1]].codClass;
                
                    id2  = allAtoms[iAN->atoms[2]].id;
                    id3  = allAtoms[iAN->atoms[1]].id;
                }
            }
            else
            {
                ha2  =allAtoms[iAN->atoms[2]].hashingValue;
                ha3  =allAtoms[iAN->atoms[1]].hashingValue;
           
                a2NB2= allAtoms[iAN->atoms[2]].codNB2Symb;
                a3NB2= allAtoms[iAN->atoms[1]].codNB2Symb;
          
                a2NB = allAtoms[iAN->atoms[2]].codNBSymb;
                a3NB = allAtoms[iAN->atoms[1]].codNBSymb;
         
                a2C  = allAtoms[iAN->atoms[2]].codClass;
                a3C  = allAtoms[iAN->atoms[1]].codClass;
                
                id2  = allAtoms[iAN->atoms[2]].id;
                id3  = allAtoms[iAN->atoms[1]].id;
                
            }
            
            std::cout << "Atom2 " <<  id2 
                      << " Its Cod class " <<  a2C << " and Hashing " << ha2 <<std::endl;
            std::cout << " Its codNBSymb " <<  a2NB << " its codNB2Symb " << a2NB2 << std::endl;
            
            std::cout << "Atom3 " <<  id3 
                      << " Its Cod class " <<  a3C << " and Hashing " << ha3 <<std::endl;
            std::cout << " Its codNBSymb " <<  a3NB << " its codNB2Symb " << a3NB2 << std::endl;
           
           
            
            int dLev = 0;
            
                    
            if((int)allDictAnglesIdx[ha1][ha2][ha3].size() !=0)
            {
                    
                std::map<ID,  std::map<ID,  std::map<ID, 
                std::map<ID,  std::map<ID,  std::map<ID,
                std::map<ID,  std::map<ID,  std::map<ID,
                int > > > > > > > > >::iterator  iFind1 
                =allDictAnglesIdx[ha1][ha2][ha3].find(a1NB2);
                
                
                if (iFind1 !=allDictAnglesIdx[ha1][ha2][ha3].end())
                {
                    std::cout << "Found " << a1NB2 << std::endl;
                }
                else
                {
                    std::cout << "not Found " << a1NB2 << " iFind1 " << std::endl;
                }
                
                
                if(iFind1 != allDictAnglesIdx[ha1][ha2][ha3].end())
                {
                    // a1NB2 matches
                    std::map<ID,  std::map<ID,  std::map<ID, 
                    std::map<ID,  std::map<ID,  std::map<ID,
                    std::map<ID,  std::map<ID,
                    int > > > > > > > >::iterator iFind2=
                    allDictAnglesIdx[ha1][ha2][ha3][a1NB2].find(a2NB2);
                
                    
                    if (iFind2 !=allDictAnglesIdx[ha1][ha2][ha3][a1NB2].end())
                    {
                        std::cout << "Found " << a2NB2 << std::endl;
                    }
                    
                    else
                    {
                        std::cout << "not Found " << a2NB2 << " iFind2 " << std::endl;
                    }
                    
                    
                    if (iFind2 !=allDictAnglesIdx[ha1][ha2][ha3][a1NB2].end())
                    {
                        //a1NB2, a2NB2 match
                        std::map<ID,  std::map<ID,  std::map<ID, 
                        std::map<ID,  std::map<ID,  std::map<ID,
                        std::map<ID,  int > > > > > > >::iterator iFind3
                        = allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2].find(a3NB2);
                        
                        if (iFind3 !=allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2].end())
                        {
                            std::cout << "Found " << a3NB2 << std::endl;
                        }
                        else
                        {
                            std::cout << "not Found " << a3NB2 << " iFind3 " << std::endl;
                        }
                        
                        if (iFind3 != allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2].end())
                        {
                            // a1NB2, a2NB2 and a3NB2 match
                            std::map<ID,  std::map<ID,  std::map<ID, 
                            std::map<ID,  std::map<ID,  std::map<ID,
                            int > > > > > >::iterator iFind4  
                            =allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2].find(a1NB);
                           
                            if (iFind4 !=allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2].end())
                            {
                                std::cout << "Found " << a1NB << std::endl;
                            }
                            else
                            {
                                std::cout << "not Found " << a1NB << " iFind4 " << std::endl;
                            }
                            
                            
                            if (iFind4 !=allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2].end())
                            {
                                // a1NB2, a2NB2, a3NB2, a1NB match
                                std::map<ID,  std::map<ID,  std::map<ID, 
                                std::map<ID,  std::map<ID,
                                int > > > > >::iterator iFind5
                                =allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB].find(a2NB);
                                
                                if (iFind5 !=allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB].end())
                                {
                                    std::cout << "Found " << a2NB << std::endl;
                                }
                                else
                                {
                                    std::cout << "not Found " << a2NB << " iFind5 " << std::endl;
                                }
                                
                                
                                if (iFind5 != allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB].end())
                                {
                                    // a1NB2, a2NB2, a3NB2, a1NB, a2NB match
                                    std::map<ID,  std::map<ID,  std::map<ID, 
                                    std::map<ID, int > > > >::iterator iFind6
                                    =allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB].find(a3NB);
                                    
                                   
                                    if (iFind6 !=allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB].end())
                                    {
                                    std::cout << "Found " << a3NB << std::endl;
                                    }
                                    else
                                    {
                                    std::cout << "not Found " << a3NB << " iFind6 " << std::endl;
                                    }
                                   
                                    
                                    if (iFind6 !=allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB].end())
                                    {
                                        // a1NB2, a2NB2, a3NB2, a1NB, a2NB, a3NB match
                                        std::map<ID, std::map<ID,  std::map<ID,
                                        int > > >::iterator iFind7 = 
                                        allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB].find(a1C);
                                       
                                        if (iFind7 !=allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB].end())
                                        {
                                            std::cout << "Found " << a1C << std::endl;
                                        }
                                        else
                                        {
                                            std::cout << "not Found " << a1C << " iFind7 "<< std::endl;
                                        }
                                        
                                        
                                        if(iFind7 != allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB].end())
                                        {
                                            // a1NB2, a2NB2, a3NB2, a1NB, a2NB, a3NB, a1C match
                                            std::map<ID, std::map<ID, int > >::iterator iFind8 =
                                            allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C].find(a2C);
                                            
                                            
                                            if (iFind8 !=allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C].end())
                                            {
                                                std::cout << "Found " << a2C << std::endl;
                                            }
                                            else
                                            {
                                                std::cout << "not Found " << a2C << " iFind8 " << std::endl;
                                            }
                                             
                                            if(iFind8 != allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C].end())
                                            {
                                                // a1NB2, a2NB2, a3NB2, a1NB, a2NB, a3NB, a1C, a2C match
                                                std::map<ID, int >::iterator iFind9 =
                                                allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C][a2C].find(a3C);
                                                
                                                
                                                if (iFind9 !=allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C][a2C].end())
                                                {
                                                std::cout << "Found " << a3C << std::endl;
                                                }
                                                else
                                                {
                                                std::cout << "not Found " << a3C << " iFind9 " << std::endl;
                                                  // exit(1);
                                                }
                                                 
                                                
                                                if(iFind9 !=
                                                allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C][a2C].end())
                                                {
                                                    // a1NB2, a2NB2, a3NB2, a1NB, a2NB, a3NB, a1C, a2C and a3C all match 
                                                    // COD has such an angle value
                                                   
                                                    int iPos =        
                                                    allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C][a2C][a3C];
                                                    if (allDictAngles[iPos].numCodValues > 5)
                                                    {
                                                        iAN->value        = allDictAngles[iPos].value;
                                                        iAN->sigValue     = allDictAngles[iPos].sigValue;
                                                        iAN->numCodValues = allDictAngles[iPos].numCodValues;
                                                        iAN->hasCodValue  = true;
                                                        std::cout << "COD finds the exact value " << iAN->value << std::endl;
                                                    }
                                                    else 
                                                    {   
                                                        std::cout << "COD has the exact matching, the account is less than 5. Using parent meaning "
                                                                  << std::endl;
                                                        std::vector<AngleDict> tDictANs;
                                                        for (std::map<ID, std::map<ID, std::map <ID, int> > >::iterator iDictANs1 =
                                                        allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB].begin();
                                                        iDictANs1 !=
                                                        allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB].end();
                                                        iDictANs1++)
                                                        {
                                                            for (std::map<ID, std::map<ID, int > >::iterator iDictANs2
                                                            =iDictANs1->second.begin();
                                                            iDictANs2 !=iDictANs1->second.end(); iDictANs2++)
                                                            {
                                                            for (std::map<ID, int >::iterator iDictANs3 
                                                            =iDictANs2->second.begin();
                                                            iDictANs3 !=iDictANs2->second.end(); iDictANs3++)
                                                            {
                                                            tDictANs.push_back(allDictAngles[iDictANs3->second]);
                                                            }
                                                            }
                                                        }
                                                        setupTargetAngleUsingMean(tDictANs, iAN);
                                                    }
                                                }
                                                else // iFind9
                                                {
                                                    // Only one outer atom COD code does not match, using the mean value
                                                    // of bonds in the group 
                                                    // allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C][a2C]
                                                    std::vector<AngleDict> tDictANs;
                                                    for (std::map<ID, int >::iterator iDictANs =
                                                         allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C][a2C].begin();
                                                         iDictANs !=
                                                         allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C][a2C].end();
                                                         iDictANs++)
                                                    {
                                                        tDictANs.push_back(allDictAngles[iDictANs->second]);
                                                    }
                                                    
                                                    setupTargetAngleUsingMean(tDictANs, iAN);
                                                }
                                            }
                                            else // iFind8
                                            {
                                                // Only two outer atom COD code do not match exactly, can still using 
                                                // the mean value of bonds in the group 
                                                // allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C]
                                                std::vector<AngleDict> tDictANs;
                                                for (std::map<ID, std::map<ID, int > >::iterator iDictANs1 =
                                                     allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C].begin();
                                                     iDictANs1 !=
                                                     allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C].end();
                                                     iDictANs1++)
                                                {
                                                    for (std::map<ID, int >::iterator iDictANs2 =iDictANs1->second.begin();
                                                         iDictANs2 !=iDictANs1->second.end(); iDictANs2++)
                                                    {
                                                        tDictANs.push_back(allDictAngles[iDictANs2->second]);
                                                    }
                                                }
                                                setupTargetAngleUsingMean(tDictANs, iAN);
                                            }
                                        }
                                        else //iFind7
                                        {
                                            // All three COD classes of component atoms do not match, still try using
                                            // the mean value of atom group
                                            // allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB]
                                            std::vector<AngleDict> tDictANs;
                                            for (std::map<ID, std::map<ID, std::map <ID, int> > >::iterator iDictANs1 =
                                                 allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB].begin();
                                                 iDictANs1 !=
                                                 allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB].end();
                                                 iDictANs1++)
                                            {
                                                for (std::map<ID, std::map<ID, int > >::iterator iDictANs2
                                                     =iDictANs1->second.begin();
                                                     iDictANs2 !=iDictANs1->second.end(); iDictANs2++)
                                                {
                                                    for (std::map<ID, int >::iterator iDictANs3 
                                                         =iDictANs2->second.begin();
                                                         iDictANs3 !=iDictANs2->second.end(); iDictANs3++)
                                                    {
                                                        tDictANs.push_back(allDictAngles[iDictANs3->second]);
                                                    }
                                                }
                                            }
                                            setupTargetAngleUsingMean(tDictANs, iAN);
                                        }
                                        
                                    }
                                    else // iFind6
                                    {
                                        // Now not only three COD classes not matching, one secondary NB atom group 
                                        // has different confs as well. Now try differently finding 
                                        // the bond with shortest atomic distances (or the lowest substitute costs)
                                        // Using the mean value could lead to weird target values because relatively
                                        // large number of samples included in calculations.
                                        std::vector<AngleDict> tDictANs;
                                        for (std::map<ID, std::map<ID, std::map <ID, std::map<ID, int> > > >::iterator iDictANs1 
                                             = allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB].begin();
                                             iDictANs1 !=
                                             allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB].end();
                                             iDictANs1++)
                                        {
                                            for(std::map<ID, std::map <ID, std::map<ID, int > > >::iterator iDictANs2
                                                = iDictANs1->second.begin(); 
                                                iDictANs2 != iDictANs1->second.end(); iDictANs2++)
                                            {
                                                for (std::map<ID, std::map<ID, int > >::iterator iDictANs3
                                                     =iDictANs2->second.begin();
                                                     iDictANs3 !=iDictANs2->second.end(); iDictANs3++)
                                                {
                                                    for (std::map<ID, int >::iterator iDictANs4
                                                         =iDictANs3->second.begin();
                                                         iDictANs4 !=iDictANs3->second.end(); iDictANs4++)
                                                    {
                                                        tDictANs.push_back(allDictAngles[iDictANs4->second]);
                                                    }
                                                }
                                            }
                                        }
                                        dLev = 1;
                                        setupTargetAngleUsingdist(tDictANs, iAN, dLev);
                                    }
                                }
                                else // iFind5
                                {
                                    // Now try differently finding he bond with shortest atomic distances
                                    // or the lowest substitute costs). Using the mean value
                                    // could lead to weird target values.
                            
                                    std::vector<AngleDict> tDictANs;
                                    for (std::map<ID, std::map<ID, std::map <ID, 
                                         std::map<ID, std::map<ID, int> > > > >::iterator iDictANs1 
                                         = allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB].begin();
                                         iDictANs1 !=
                                         allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB].end();
                                         iDictANs1++)
                                    {
                                        for(std::map<ID, std::map <ID, 
                                            std::map<ID, std::map<ID, int > > > >::iterator iDictANs2
                                            = iDictANs1->second.begin(); 
                                            iDictANs2 != iDictANs1->second.end(); iDictANs2++)
                                        {
                                            for (std::map<ID, std::map<ID, std::map<ID, int > > >::iterator iDictANs3
                                                 =iDictANs2->second.begin();
                                                 iDictANs3 !=iDictANs2->second.end(); iDictANs3++)
                                            {
                                                for (std::map<ID, std::map<ID, int > >::iterator iDictANs4
                                                     =iDictANs3->second.begin();
                                                     iDictANs4 !=iDictANs3->second.end(); iDictANs4++)
                                                {
                                                    for (std::map<ID, int >::iterator iDictANs5
                                                         =iDictANs4->second.begin();
                                                         iDictANs5 !=iDictANs4->second.end(); iDictANs5++)
                                                    {
                                                        tDictANs.push_back(allDictAngles[iDictANs5->second]);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    // std::cout << "exit(1) here " << std::endl;
                                    dLev = 1;
                                    setupTargetAngleUsingdist(tDictANs, iAN, dLev);    
                                }
                            }
                            else // iFind4
                            {
                                // Now try differently finding the bond with shortest atomic distances
                                // or the lowest substitute costs). Using the mean value
                                // could lead to weird target values.
                                
                                std::vector<AngleDict> tDictANs;
                                for (std::map<ID, std::map<ID, std::map <ID, std::map<ID, 
                                     std::map<ID, std::map<ID, int> > > > > >::iterator iDictANs1 
                                     = allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2].begin();
                                     iDictANs1 !=
                                     allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2].end();
                                     iDictANs1++)
                                {
                                    for(std::map<ID, std::map <ID, std::map<ID,  
                                        std::map<ID, std::map<ID, int > > > > >::iterator iDictANs2
                                        = iDictANs1->second.begin(); 
                                       iDictANs2 != iDictANs1->second.end(); iDictANs2++)
                                    {
                                        for(std::map<ID, std::map <ID, 
                                            std::map<ID, std::map<ID, int > > > >::iterator iDictANs3
                                            = iDictANs2->second.begin(); 
                                            iDictANs3 != iDictANs2->second.end(); iDictANs3++)
                                        {
                                            for (std::map<ID, std::map<ID, std::map<ID, int> > >::iterator iDictANs4
                                                 =iDictANs3->second.begin();
                                                 iDictANs4 !=iDictANs3->second.end(); iDictANs4++)
                                            {
                                                for (std::map<ID, std::map<ID, int > >::iterator iDictANs5
                                                     =iDictANs4->second.begin();
                                                     iDictANs5 !=iDictANs4->second.end(); iDictANs5++)
                                                {
                                                    for (std::map<ID, int >::iterator iDictANs6
                                                         =iDictANs5->second.begin();
                                                         iDictANs6 !=iDictANs5->second.end(); iDictANs6++)
                                                    {
                                                        std::cout << "angle idx " << iDictANs6->second << std::endl;
                                                        std::cout << "angle value  " << (int)allDictAngles[iDictANs6->second].value 
                                                                  << std::endl; 
                                                        tDictANs.push_back(allDictAngles[iDictANs6->second]);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                
                                dLev = 1;
                                setupTargetAngleUsingdist(tDictANs, iAN, dLev);   
                            }
                        }
                        else // iFind3
                        {
                            // Finding the bond with shortest atomic distances 
                            // (or the lowest substitute costs)
                            std::vector<AngleDict> tDictANs;
                            for (std::map<ID, std::map<ID, std::map <ID, std::map<ID, std::map<ID, 
                                 std::map<ID, std::map<ID, int> > > > > > >::iterator iDictANs1 
                                 = allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2].begin();
                                 iDictANs1 !=
                                 allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2].end();
                                 iDictANs1++)
                            {
                                for(std::map<ID, std::map <ID, std::map<ID, std::map<ID, 
                                    std::map<ID, std::map<ID, int > > > > > >::iterator iDictANs2
                                    = iDictANs1->second.begin(); 
                                    iDictANs2 != iDictANs1->second.end(); iDictANs2++)
                                {
                                    for(std::map<ID, std::map <ID, std::map<ID,  
                                        std::map<ID, std::map<ID, int > > > > >::iterator iDictANs3
                                        = iDictANs2->second.begin(); 
                                       iDictANs3 != iDictANs2->second.end(); iDictANs3++)
                                    {
                                        for(std::map<ID, std::map <ID, 
                                            std::map<ID, std::map<ID, int > > > >::iterator iDictANs4
                                            = iDictANs3->second.begin(); 
                                            iDictANs4 != iDictANs3->second.end(); iDictANs4++)
                                        {
                                            for (std::map<ID, std::map<ID, std::map<ID, int > > >::iterator iDictANs5
                                                 =iDictANs4->second.begin();
                                                 iDictANs5 !=iDictANs4->second.end(); iDictANs5++)
                                            {
                                                for (std::map<ID, std::map<ID, int > >::iterator iDictANs6
                                                     =iDictANs5->second.begin();
                                                     iDictANs6 !=iDictANs5->second.end(); iDictANs6++)
                                                {
                                                    for (std::map<ID, int >::iterator iDictANs7
                                                         =iDictANs6->second.begin();
                                                         iDictANs7 !=iDictANs6->second.end(); iDictANs7++)
                                                    {
                                                        tDictANs.push_back(allDictAngles[iDictANs7->second]);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            dLev = 2;
                            setupTargetAngleUsingdist(tDictANs, iAN, dLev);
                        }
                    }
                    else // iFind2 
                    {
                        // Finding the bond with shortest atomic distances 
                        // (or the lowest substitute costs)
                        std::vector<AngleDict> tDictANs;
                        for (std::map<ID, std::map<ID, std::map <ID, std::map<ID,  
                             std::map<ID, std::map<ID, std::map<ID, 
                             std::map<ID, int> > > > > > > >::iterator iDictANs1 
                             = allDictAnglesIdx[ha1][ha2][ha3][a1NB2].begin();
                             iDictANs1 !=
                             allDictAnglesIdx[ha1][ha2][ha3][a1NB2].end();
                                 iDictANs1++)
                        {
                            for(std::map<ID, std::map <ID, std::map<ID, 
                                std::map<ID, std::map<ID, std::map<ID,  
                                std::map<ID, int > > > > > > >::iterator iDictANs2
                                    = iDictANs1->second.begin(); 
                                    iDictANs2 != iDictANs1->second.end(); iDictANs2++)
                            {
                                for(std::map<ID, std::map <ID, std::map<ID, std::map<ID, 
                                    std::map<ID, std::map<ID, int > > > > > >::iterator iDictANs3
                                        = iDictANs2->second.begin(); 
                                    iDictANs3 != iDictANs2->second.end(); iDictANs3++)
                                {
                                    for(std::map<ID, std::map <ID, std::map<ID,  
                                        std::map<ID, std::map<ID, int > > > > >::iterator iDictANs4
                                        = iDictANs3->second.begin(); 
                                       iDictANs4 != iDictANs3->second.end(); iDictANs4++)
                                    {
                                        for(std::map<ID, std::map <ID, 
                                            std::map<ID, std::map<ID, int > > > >::iterator iDictANs5
                                            = iDictANs4->second.begin(); 
                                            iDictANs5 != iDictANs4->second.end(); iDictANs5++)
                                        {
                                            for (std::map<ID, std::map<ID, std::map<ID, int > > >::iterator iDictANs6
                                                 =iDictANs5->second.begin();
                                                 iDictANs6 !=iDictANs5->second.end(); iDictANs6++)
                                            {
                                                for (std::map<ID, std::map<ID, int > >::iterator iDictANs7
                                                     =iDictANs6->second.begin();
                                                     iDictANs7 !=iDictANs6->second.end(); iDictANs7++)
                                                {
                                                    for (std::map<ID, int >::iterator iDictANs8
                                                         =iDictANs7->second.begin();
                                                         iDictANs8 !=iDictANs7->second.end(); iDictANs8++)
                                                    {
                                                        tDictANs.push_back(allDictAngles[iDictANs8->second]);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        dLev = 2;
                        setupTargetAngleUsingdist(tDictANs, iAN, dLev);
                    }
                }
                else // iFind1
                {
                    // Finding the bond with shortest atomic distances 
                    // (or the lowest substitute costs)
                    std::vector<AngleDict> tDictANs;
                    for (std::map<ID, std::map<ID, std::map <ID, std::map<ID,  
                         std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                         std::map<ID, int> > > > > > > > >::iterator iDictANs1 
                             = allDictAnglesIdx[ha1][ha2][ha3].begin();
                         iDictANs1 !=
                             allDictAnglesIdx[ha1][ha2][ha3].end();
                         iDictANs1++)
                    {
                        for(std::map<ID, std::map <ID, std::map<ID, 
                            std::map<ID, std::map<ID, std::map<ID, std::map<ID,  
                            std::map<ID, int > > > > > > > >::iterator iDictANs2
                            = iDictANs1->second.begin(); 
                            iDictANs2 != iDictANs1->second.end(); iDictANs2++)
                        {
                            for(std::map<ID, std::map <ID, std::map<ID, 
                                std::map<ID, std::map<ID, std::map<ID,  
                                std::map<ID, int > > > > > > >::iterator iDictANs3
                                    = iDictANs2->second.begin(); 
                                iDictANs3 != iDictANs2->second.end(); iDictANs3++)
                            {
                                for(std::map<ID, std::map <ID, std::map<ID, std::map<ID, 
                                    std::map<ID, std::map<ID, int > > > > > >::iterator iDictANs4
                                        = iDictANs3->second.begin(); 
                                    iDictANs4 != iDictANs3->second.end(); iDictANs4++)
                                {
                                    for(std::map<ID, std::map <ID, std::map<ID,  
                                        std::map<ID, std::map<ID, int > > > > >::iterator iDictANs5
                                        = iDictANs4->second.begin(); 
                                       iDictANs5 != iDictANs4->second.end(); iDictANs5++)
                                    {
                                        for(std::map<ID, std::map <ID, 
                                            std::map<ID, std::map<ID, int > > > >::iterator iDictANs6
                                            = iDictANs5->second.begin(); 
                                            iDictANs6 != iDictANs5->second.end(); iDictANs6++)
                                        {
                                            for (std::map<ID, std::map<ID, 
                                                 std::map<ID, int > > >::iterator iDictANs7
                                                 =iDictANs6->second.begin();
                                                 iDictANs7 !=iDictANs6->second.end(); iDictANs7++)
                                            {
                                                for (std::map<ID, std::map<ID, int > >::iterator iDictANs8
                                                     =iDictANs7->second.begin();
                                                     iDictANs8 !=iDictANs7->second.end(); iDictANs8++)
                                                {
                                                    for (std::map<ID, int >::iterator iDictANs9
                                                         =iDictANs8->second.begin();
                                                         iDictANs9 !=iDictANs8->second.end(); iDictANs9++)
                                                    {
                                                        tDictANs.push_back(allDictAngles[iDictANs9->second]);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    dLev = 2;
                    setupTargetAngleUsingdist(tDictANs, iAN, dLev);
                }
                //std::cout << "Final angle is " << iAN->value << std::endl;
                //if (iAN->hasCodValue)
                //{
                //    std::cout << "It has exact matches to a COD angle " << std::endl;
                //}
                //else
                //{
                //    std::cout << "It is generated by either averaging or shortest distance " << std::endl;
                //}
            }
            else   
            {
                // could not find three exact matches on 3 atomic hashing values
                // using approximate default values
                if (allAtoms[iAN->atoms[0]].bondingIdx >0 && 
                    allAtoms[iAN->atoms[0]].bondingIdx <4)
                {
                    // std::cout << "Center atom bond index is  " << allAtoms[iAN->atoms[0]].bondingIdx<< std::endl;
                    iAN->value = DefaultOrgAngles[allAtoms[iAN->atoms[0]].bondingIdx];
                }
                else
                {
                    std::cout << "Could not find COD angle value and even default angle value."
                        <<std::endl << "The atoms in the target angle are: " << std::endl
                        << "The inner atom "   << allAtoms[iAN->atoms[0]].id << std::endl
                        << "Its COD class is " << allAtoms[iAN->atoms[0]].codClass << std::endl
                        << "The outer atom 1 " << allAtoms[iAN->atoms[1]].id << std::endl
                        << "Its COD class is " << allAtoms[iAN->atoms[1]].codClass << std::endl
                        << "The outer atom 2 "  << allAtoms[iAN->atoms[2]].id << std::endl
                        << "Its COD class is " << allAtoms[iAN->atoms[2]].codClass << std::endl;
                
                    exit(1);
                }
                 
            }
        
    }
    
    void CodClassify::searchCodOrgAngles2(std::vector<AngleDict>::iterator iAN)
    {    
            /*
            codClassToAtomAng(allAtoms[iAN->atoms[0]].codClass, allAtoms[iAN->atoms[0]]);
            codClassToAtomAng(allAtoms[iAN->atoms[1]].codClass, allAtoms[iAN->atoms[1]]);
            codClassToAtomAng(allAtoms[iAN->atoms[2]].codClass, allAtoms[iAN->atoms[2]]);
            */
        
            int ha1, ha2, ha3, s2, s3;
            ID a1NB2, a1NB, a1C, a2NB2, a2NB, a2C,a3NB2, a3NB, a3C;
            ID id2, id3;
            ha1   = allAtoms[iAN->atoms[0]].hashingValue;
            a1NB2 = allAtoms[iAN->atoms[0]].codNB2Symb;
            a1NB  = allAtoms[iAN->atoms[0]].codNBSymb;
            a1C   = allAtoms[iAN->atoms[0]].codClass;
            
            std::cout << "Atom1 " <<  allAtoms[iAN->atoms[0]].id 
                      << " Its Cod class " <<  a1C <<  " and Hashing "
                      << ha1 <<std::endl;
            std::cout << " Its codNBSymb " <<  a1NB << " its codNB2Symb " << a1NB2 << std::endl;
            
            if ((int)allAtoms[iAN->atoms[1]].hashingValue
                    <(int)allAtoms[iAN->atoms[2]].hashingValue )
            {
                s2 =1;
                s3 =2;
                ha2  =allAtoms[iAN->atoms[1]].hashingValue;
                ha3  =allAtoms[iAN->atoms[2]].hashingValue;
           
                a2NB2= allAtoms[iAN->atoms[1]].codNB2Symb;
                a3NB2= allAtoms[iAN->atoms[2]].codNB2Symb;
          
                a2NB = allAtoms[iAN->atoms[1]].codNBSymb;
                a3NB = allAtoms[iAN->atoms[2]].codNBSymb;
         
                a2C  = allAtoms[iAN->atoms[1]].codClass;
                a3C  = allAtoms[iAN->atoms[2]].codClass;
                
                id2  = allAtoms[iAN->atoms[1]].id;
                id3  = allAtoms[iAN->atoms[2]].id;
            }
            else if ((int)allAtoms[iAN->atoms[1]].hashingValue
                    ==(int)allAtoms[iAN->atoms[2]].hashingValue)
            {
                if ((int)allAtoms[iAN->atoms[1]].codClass.size() <=
                        (int)allAtoms[iAN->atoms[2]].codClass.size())
                {
                    s2 =1;
                    s3 =2;
                    ha2  =allAtoms[iAN->atoms[1]].hashingValue;
                    ha3  =allAtoms[iAN->atoms[2]].hashingValue;
           
                    a2NB2= allAtoms[iAN->atoms[1]].codNB2Symb;
                    a3NB2= allAtoms[iAN->atoms[2]].codNB2Symb;
          
                    a2NB = allAtoms[iAN->atoms[1]].codNBSymb;
                    a3NB = allAtoms[iAN->atoms[2]].codNBSymb;
         
                    a2C  = allAtoms[iAN->atoms[1]].codClass;
                    a3C  = allAtoms[iAN->atoms[2]].codClass;
                
                    id2  = allAtoms[iAN->atoms[1]].id;
                    id3  = allAtoms[iAN->atoms[2]].id;
                }
                else
                {
                    s2 =2;
                    s3 =1;
                    ha2  =allAtoms[iAN->atoms[2]].hashingValue;
                    ha3  =allAtoms[iAN->atoms[1]].hashingValue;
           
                    a2NB2= allAtoms[iAN->atoms[2]].codNB2Symb;
                    a3NB2= allAtoms[iAN->atoms[1]].codNB2Symb;
          
                    a2NB = allAtoms[iAN->atoms[2]].codNBSymb;
                    a3NB = allAtoms[iAN->atoms[1]].codNBSymb;
         
                    a2C  = allAtoms[iAN->atoms[2]].codClass;
                    a3C  = allAtoms[iAN->atoms[1]].codClass;
                
                    id2  = allAtoms[iAN->atoms[2]].id;
                    id3  = allAtoms[iAN->atoms[1]].id;
                }
            }
            else
            {
                s2 =2;
                s3 =1;
                ha2  =allAtoms[iAN->atoms[2]].hashingValue;
                ha3  =allAtoms[iAN->atoms[1]].hashingValue;
           
                a2NB2= allAtoms[iAN->atoms[2]].codNB2Symb;
                a3NB2= allAtoms[iAN->atoms[1]].codNB2Symb;
          
                a2NB = allAtoms[iAN->atoms[2]].codNBSymb;
                a3NB = allAtoms[iAN->atoms[1]].codNBSymb;
         
                a2C  = allAtoms[iAN->atoms[2]].codClass;
                a3C  = allAtoms[iAN->atoms[1]].codClass;
                
                id2  = allAtoms[iAN->atoms[2]].id;
                id3  = allAtoms[iAN->atoms[1]].id;
                
            }
            
            std::cout << "Atom2 " <<  id2 
                      << " Its Cod class " <<  a2C << " and Hashing " << ha2 <<std::endl;
            std::cout << " Its codNBSymb " <<  a2NB << " its codNB2Symb " << a2NB2 << std::endl;
            
            std::cout << "Atom3 " <<  id3 
                      << " Its Cod class " <<  a3C << " and Hashing " << ha3 <<std::endl;
            std::cout << " Its codNBSymb " <<  a3NB << " its codNB2Symb " << a3NB2 << std::endl;
           
           
            
            int dLev = 0;
            
            
          
                    
            if(allDictAnglesIdx2.find(ha1) != allDictAnglesIdx2.end() &&
               allDictAnglesIdx2[ha1].find(ha2) != allDictAnglesIdx2[ha1].end() &&
               allDictAnglesIdx2[ha1][ha2].find(ha3) != allDictAnglesIdx2[ha1][ha2].end())
            {
                    
                std::map<ID,  std::map<ID,  std::map<ID, 
                std::vector<aValueSet> > > >::iterator  iFind1 
                =allDictAnglesIdx2[ha1][ha2][ha3].find(a1NB2);
                
                
                if (iFind1 !=allDictAnglesIdx2[ha1][ha2][ha3].end())
                {
                    std::cout << "Found " << a1NB2 << std::endl;
                }
                else
                {
                    std::cout << "not Found " << a1NB2 << " iFind1 " << std::endl;
                }
               
                
                if(iFind1 != allDictAnglesIdx2[ha1][ha2][ha3].end())
                {
                    // a1NB2 matches
                    std::map<ID,  std::map<ID,
                    std::vector<aValueSet> > >::iterator iFind2=
                    allDictAnglesIdx2[ha1][ha2][ha3][a1NB2].find(a2NB2);
                
                     
                    if (iFind2 !=allDictAnglesIdx2[ha1][ha2][ha3][a1NB2].end())
                    {
                        std::cout << "Found " << a2NB2 << std::endl;
                    }
                    
                    else
                    {
                        std::cout << "not Found " << a2NB2 << " iFind2 " << std::endl;
                    }
                    
                   
                    
                    if (iFind2 !=allDictAnglesIdx2[ha1][ha2][ha3][a1NB2].end())
                    {
                        //a1NB2, a2NB2 match
                       
                        std::map<ID,  std::vector<aValueSet> >::iterator iFind3
                        = allDictAnglesIdx2[ha1][ha2][ha3][a1NB2][a2NB2].find(a3NB2);
                        
                        
                        if (iFind3 !=allDictAnglesIdx2[ha1][ha2][ha3][a1NB2][a2NB2].end())
                        {
                            std::cout << "Found " << a3NB2 << std::endl;
                        }
                        else
                        {
                            std::cout << "not Found " << a3NB2 << " iFind3 " << std::endl;
                        }
                        
                        
                        if (iFind3 != allDictAnglesIdx2[ha1][ha2][ha3][a1NB2][a2NB2].end())
                        {
                            // a1NB2, a2NB2 and a3NB2 match
                            std::map<ID,  std::map<ID,  std::map<ID, 
                            std::map<ID,  std::map<ID,  std::map<ID,
                            int > > > > > >::iterator iFind4  
                            =allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2].find(a1NB);
                            
                            
                            if (iFind4 !=allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2].end())
                            {
                                std::cout << "Found " << a1NB << std::endl;
                            }
                            else
                            {
                                std::cout << "not Found " << a1NB << " iFind4 " << std::endl;
                            }
                            
                            
                            if (iFind4 !=allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2].end())
                            {
     
                                // a1NB2, a2NB2, a3NB2, a1NB match
                                std::map<ID,  std::map<ID,  std::map<ID, 
                                std::map<ID,  std::map<ID,
                                int > > > > >::iterator iFind5
                                =allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB].find(a2NB);
                                
                                
                                if (iFind5 !=allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB].end())
                                {
                                    std::cout << "Found " << a2NB << std::endl;
                                }
                                else
                                {
                                    std::cout << "not Found " << a2NB << " iFind5 " << std::endl;
                                }
                                
                                
                                if (iFind5 != allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB].end())
                                {
                                    // a1NB2, a2NB2, a3NB2, a1NB, a2NB match
                                    std::map<ID,  std::map<ID,  std::map<ID, 
                                    std::map<ID, int > > > >::iterator iFind6
                                    =allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB].find(a3NB);
                                    
                                    
                                    if (iFind6 !=allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB].end())
                                    {
                                    std::cout << "Found " << a3NB << std::endl;
                                    }
                                    else
                                    {
                                    std::cout << "not Found " << a3NB << " iFind6 " << std::endl;
                                    }
                                    
                                    
                                    if (iFind6 !=allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB].end())
                                    {
                                        // a1NB2, a2NB2, a3NB2, a1NB, a2NB, a3NB match
                                        std::map<ID, std::map<ID,  std::map<ID,
                                        int > > >::iterator iFind7 = 
                                        allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB].find(a1C);
                                        
                                        
                                        if (iFind7 !=allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB].end())
                                        {
                                            std::cout << "Found " << a1C << std::endl;
                                        }
                                        else
                                        {
                                            std::cout << "not Found " << a1C << " iFind7 "<< std::endl;
                                        }
                                        
                                        
                                        if(iFind7 != allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB].end())
                                        {
                                            // a1NB2, a2NB2, a3NB2, a1NB, a2NB, a3NB, a1C match
                                            std::map<ID, std::map<ID, int > >::iterator iFind8 =
                                            allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C].find(a2C);
                                            
                                            
                                            if (iFind8 !=allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C].end())
                                            {
                                                std::cout << "Found " << a2C << std::endl;
                                            }
                                            else
                                            {
                                                std::cout << "not Found " << a2C << " iFind8 " << std::endl;
                                            }
                                            
                                             
                                            if(iFind8 != allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C].end())
                                            {
                                                // a1NB2, a2NB2, a3NB2, a1NB, a2NB, a3NB, a1C, a2C match
                                                std::map<ID, int >::iterator iFind9 =
                                                allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C][a2C].find(a3C);
                                                
                                                /*
                                                if (iFind9 !=allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C][a2C].end())
                                                { 
                                                   std::cout << "Found " << a3C << std::endl;
                                                }
                                                else
                                                {
                                                    std::cout << "not Found " << a3C << " iFind9 " << std::endl;
                                                  // exit(1);
                                                }
                                                */
                                                
                                                if(iFind9 !=
                                                allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C][a2C].end())
                                                {
                                                    // a1NB2, a2NB2, a3NB2, a1NB, a2NB, a3NB, a1C, a2C and a3C all match 
                                                    // COD has such an angle value
                                                          
                                                    int iPos =        
                                                    allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C][a2C][a3C];
                                                    //std::cout << "iPos =" << iPos << std::endl;
                                                    iAN->value        = allDictAngles[iPos].value;
                                                    iAN->sigValue     = allDictAngles[iPos].sigValue;
                                                    iAN->numCodValues = allDictAngles[iPos].numCodValues;
                                                    iAN->hasCodValue  = true;
                                                    iAN->levelCodValue = 0;
                                                    std::cout << "COD finds the exact value " << iAN->value << std::endl;
                                                    std::cout << "the sigma value " << iAN->sigValue << std::endl;
                                                    
                                                }
                                                else // iFind9
                                                {
                                                    /*
                                                    std::vector<AngleDict> tDictANs;
                                                    for (std::map<ID, int >::iterator iDictANs =
                                                         allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C][a2C].begin();
                                                         iDictANs !=
                                                         allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C][a2C].end();
                                                         iDictANs++)
                                                    {
                                                        tDictANs.push_back(allDictAngles[iDictANs->second]);
                                                    }
                                                    
                                                    setupTargetAngleUsingMean2(tDictANs, iAN, ha1, ha2, ha3, a1NB2, a2NB2, a3NB2);
                                                    */
                                                    iAN->value        = 
                                                    allDictAnglesIdx1[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][0].value;
                                                    iAN->sigValue     =
                                                    allDictAnglesIdx1[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][0].sigValue;
                                                    iAN->numCodValues = 
                                                    allDictAnglesIdx1[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][0].numCodValues;
                                                    iAN->levelCodValue = 1;
                                                }
                                            }
                                            else // iFind8
                                            {
                                                // Only two outer atom COD code do not match exactly, can still using 
                                                // the mean value of bonds in the group 
                                                // allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C]
                                                /*
                                                std::vector<AngleDict> tDictANs;
                                                for (std::map<ID, std::map<ID, int > >::iterator iDictANs1 =
                                                     allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C].begin();
                                                     iDictANs1 !=
                                                     allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1C].end();
                                                     iDictANs1++)
                                                {
                                                    for (std::map<ID, int >::iterator iDictANs2 =iDictANs1->second.begin();
                                                         iDictANs2 !=iDictANs1->second.end(); iDictANs2++)
                                                    {
                                                        tDictANs.push_back(allDictAngles[iDictANs2->second]);
                                                    }
                                                 }
                                                 
                                                setupTargetAngleUsingMean2(tDictANs, iAN, ha1, ha2, ha3, a1NB2, a2NB2, a3NB2);
                                                */
                                                iAN->value        = 
                                                    allDictAnglesIdx1[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][0].value;
                                                iAN->sigValue     =
                                                    allDictAnglesIdx1[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][0].sigValue;
                                                iAN->numCodValues = 
                                                    allDictAnglesIdx1[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][0].numCodValues;
                                                iAN->levelCodValue = 1;
                                            }
                                        }
                                        else //iFind7
                                        {
                                            // All three COD classes of component atoms do not match, still try using
                                            // the mean value of atom group
                                            // allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB]
                                            /*
                                            std::vector<AngleDict> tDictANs;
                                            for (std::map<ID, std::map<ID, std::map <ID, int> > >::iterator iDictANs1 =
                                                 allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB].begin();
                                                 iDictANs1 !=
                                                 allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB].end();
                                                 iDictANs1++)
                                            {
                                                for (std::map<ID, std::map<ID, int > >::iterator iDictANs2
                                                     =iDictANs1->second.begin();
                                                     iDictANs2 !=iDictANs1->second.end(); iDictANs2++)
                                                {
                                                    for (std::map<ID, int >::iterator iDictANs3 
                                                         =iDictANs2->second.begin();
                                                         iDictANs3 !=iDictANs2->second.end(); iDictANs3++)
                                                    {
                                                        tDictANs.push_back(allDictAngles[iDictANs3->second]);
                                                    }
                                                }
                                            }
                                            
                                            setupTargetAngleUsingMean2(tDictANs, iAN, ha1, ha2, ha3, a1NB2, a2NB2, a3NB2);
                                            */
                                            iAN->value        = 
                                                    allDictAnglesIdx1[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][0].value;
                                            iAN->sigValue     =
                                                    allDictAnglesIdx1[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][0].sigValue;
                                            iAN->numCodValues = 
                                                    allDictAnglesIdx1[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][0].numCodValues;
                                            iAN->levelCodValue =1;
                                        }
                                        
                                    }
                                    else // iFind6
                                    {
                                        
                                        std::vector<AngleDict> tDictANs;
                                        for (std::map<ID, std::map<ID, std::map <ID, std::map<ID, int> > > >::iterator iDictANs1 
                                             = allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB].begin();
                                             iDictANs1 !=
                                             allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB].end();
                                             iDictANs1++)
                                        {
                                            for(std::map<ID, std::map <ID, std::map<ID, int > > >::iterator iDictANs2
                                                = iDictANs1->second.begin(); 
                                                iDictANs2 != iDictANs1->second.end(); iDictANs2++)
                                            {
                                                for (std::map<ID, std::map<ID, int > >::iterator iDictANs3
                                                     =iDictANs2->second.begin();
                                                     iDictANs3 !=iDictANs2->second.end(); iDictANs3++)
                                                {
                                                    for (std::map<ID, int >::iterator iDictANs4
                                                         =iDictANs3->second.begin();
                                                         iDictANs4 !=iDictANs3->second.end(); iDictANs4++)
                                                    {
                                                        tDictANs.push_back(allDictAngles[iDictANs4->second]);
                                                    }
                                                }
                                            }
                                        }

                                        setupTargetAngleUsingMean2(tDictANs, iAN, ha1, ha2, ha3, a1NB2, a2NB2, a3NB2);
                                        iAN->levelCodValue = 2;
                                        /*                                        
                                        iAN->value        = 
                                                    allDictAnglesIdx2[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][0].value;
                                        iAN->sigValue     =
                                                    allDictAnglesIdx2[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][0].sigValue;
                                        iAN->numCodValues = 
                                                    allDictAnglesIdx2[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][0].numCodValues;
                                        */
                                    }
                                }
                                else // iFind5
                                {
                                    // Now try differently finding the angle with shortest atomic distances
                                    // or the lowest substitute costs). Using the mean value
                                    // could lead to weird target values.
                                    std::vector<AngleDict> tDictANs;
                                    for (std::map<ID, std::map<ID, std::map <ID, 
                                         std::map<ID, std::map<ID, int> > > > >::iterator iDictANs1 
                                         = allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB].begin();
                                         iDictANs1 !=
                                         allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB].end();
                                         iDictANs1++)
                                    {
                                        for(std::map<ID, std::map <ID, 
                                            std::map<ID, std::map<ID, int > > > >::iterator iDictANs2
                                            = iDictANs1->second.begin(); 
                                            iDictANs2 != iDictANs1->second.end(); iDictANs2++)
                                        {
                                            for (std::map<ID, std::map<ID, std::map<ID, int > > >::iterator iDictANs3
                                                 =iDictANs2->second.begin();
                                                 iDictANs3 !=iDictANs2->second.end(); iDictANs3++)
                                            {
                                                for (std::map<ID, std::map<ID, int > >::iterator iDictANs4
                                                     =iDictANs3->second.begin();
                                                     iDictANs4 !=iDictANs3->second.end(); iDictANs4++)
                                                {
                                                    for (std::map<ID, int >::iterator iDictANs5
                                                         =iDictANs4->second.begin();
                                                         iDictANs5 !=iDictANs4->second.end(); iDictANs5++)
                                                    {
                                                        tDictANs.push_back(allDictAngles[iDictANs5->second]);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    
                                    setupTargetAngleUsingMean2(tDictANs, iAN, ha1, ha2, ha3, a1NB2, a2NB2, a3NB2);
                                    iAN->levelCodValue = 2;
                                    
                                    //iAN->value        = 
                                    //                allDictAnglesIdx2[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][0].value;
                                    //iAN->sigValue     =
                                    //                allDictAnglesIdx2[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][0].sigValue;
                                    //iAN->numCodValues = 
                                    //                allDictAnglesIdx2[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][0].numCodValues;
                                       
                                }
                            }
                            else // iFind4
                            {
                                
                                
                                
                                iAN->value        = 
                                                    allDictAnglesIdx2[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][0].value;
                                iAN->sigValue     =
                                                    allDictAnglesIdx2[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][0].sigValue;
                                iAN->numCodValues = 
                                                    allDictAnglesIdx2[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][0].numCodValues;
                                iAN->levelCodValue = 2;
                                
                            }
                        }
                        else // iFind3
                        {
                            
                            if (allDictAnglesIdx3[ha1][ha2][ha3][0].numCodValues >=10 
                                    && allDictAnglesIdx3[ha1][ha2][ha3][0].sigValue < 5.00)
                            {
                                iAN->value        = 
                                                    allDictAnglesIdx3[ha1][ha2][ha3][0].value;
                                iAN->sigValue     =
                                                    allDictAnglesIdx3[ha1][ha2][ha3][0].sigValue;
                                iAN->numCodValues = 
                                                    allDictAnglesIdx3[ha1][ha2][ha3][0].numCodValues;
                            }
                            else
                            {
                                bool aCCP4S=getCCP4Angle(iAN);
                            
                                if (!aCCP4S)
                                {
                                    // Finding the bond with shortest atomic distances 
                                    // (or the lowest substitute costs)
                            
                                    // std::vector<AngleDict> tDictANs;
                                    std::vector<REAL> tAllValues;
                                    std::vector<int>  tAllNums;
                                    for (std::map<ID, std::map<ID, std::map <ID, std::map<ID, std::map<ID, 
                                         std::map<ID, std::map<ID, int> > > > > > >::iterator iDictANs1 
                                         = allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2].begin();
                                         iDictANs1 !=
                                         allDictAnglesIdx[ha1][ha2][ha3][a1NB2][a2NB2].end();
                                         iDictANs1++)
                                    {
                                        for(std::map<ID, std::map <ID, std::map<ID, std::map<ID, 
                                            std::map<ID, std::map<ID, int > > > > > >::iterator iDictANs2
                                            = iDictANs1->second.begin(); 
                                            iDictANs2 != iDictANs1->second.end(); iDictANs2++)
                                        {
                                            for(std::map<ID, std::map <ID, std::map<ID,  
                                                std::map<ID, std::map<ID, int > > > > >::iterator iDictANs3
                                                = iDictANs2->second.begin(); 
                                                iDictANs3 != iDictANs2->second.end(); iDictANs3++)
                                            {
                                                for(std::map<ID, std::map <ID, 
                                                    std::map<ID, std::map<ID, int > > > >::iterator iDictANs4
                                                    = iDictANs3->second.begin(); 
                                                    iDictANs4 != iDictANs3->second.end(); iDictANs4++)
                                                {
                                                    for (std::map<ID, std::map<ID, std::map<ID, int > > >::iterator iDictANs5
                                                         =iDictANs4->second.begin();
                                                         iDictANs5 !=iDictANs4->second.end(); iDictANs5++)
                                                    {
                                                        for (std::map<ID, std::map<ID, int > >::iterator iDictANs6
                                                             =iDictANs5->second.begin();
                                                             iDictANs6 !=iDictANs5->second.end(); iDictANs6++)
                                                        {
                                                            for (std::map<ID, int >::iterator iDictANs7
                                                                 =iDictANs6->second.begin();
                                                                 iDictANs7 !=iDictANs6->second.end(); iDictANs7++)
                                                            {
                                                                // tDictANs.push_back(allDictAngles[iDictANs7->second]);
                                                                tAllValues.push_back(allDictAngles[iDictANs7->second].value);
                                                                tAllNums.push_back(allDictAngles[iDictANs7->second].numCodValues);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                
                                    //dLev = 3;
                                    //setupTargetAngleUsingdist2(tDictANs, iAN, s2, s3, dLev);
                                    setupTargetAngleUsingMean3(tAllValues, tAllNums, iAN);
                                }
                            }
                            iAN->levelCodValue = 3;
                        }
                    }
                    else // iFind2 
                    {
                        if (allDictAnglesIdx3[ha1][ha2][ha3][0].numCodValues >=10 
                                    && allDictAnglesIdx3[ha1][ha2][ha3][0].sigValue < 5.00)
                        {
                            iAN->value        = allDictAnglesIdx3[ha1][ha2][ha3][0].value;
                            iAN->sigValue     = allDictAnglesIdx3[ha1][ha2][ha3][0].sigValue;
                            iAN->numCodValues = allDictAnglesIdx3[ha1][ha2][ha3][0].numCodValues;
                        }
                        else
                        {
                            bool aCCP4S=getCCP4Angle(iAN);
                            
                            if (!aCCP4S)
                            {
                                // Finding the bond with shortest atomic distances 
                                // (or the lowest substitute costs)
                                std::vector<AngleDict> tDictANs;
                                for (std::map<ID, std::map<ID, std::map <ID, std::map<ID,  
                                     std::map<ID, std::map<ID, std::map<ID, 
                                     std::map<ID, int> > > > > > > >::iterator iDictANs1 
                                     = allDictAnglesIdx[ha1][ha2][ha3][a1NB2].begin();
                                     iDictANs1 !=
                                     allDictAnglesIdx[ha1][ha2][ha3][a1NB2].end();
                                     iDictANs1++)
                                {
                                    for(std::map<ID, std::map <ID, std::map<ID, 
                                        std::map<ID, std::map<ID, std::map<ID,  
                                        std::map<ID, int > > > > > > >::iterator iDictANs2
                                        = iDictANs1->second.begin(); 
                                        iDictANs2 != iDictANs1->second.end(); iDictANs2++)
                                    {
                                        for(std::map<ID, std::map <ID, std::map<ID, std::map<ID, 
                                            std::map<ID, std::map<ID, int > > > > > >::iterator iDictANs3
                                            = iDictANs2->second.begin(); 
                                            iDictANs3 != iDictANs2->second.end(); iDictANs3++)
                                        {
                                            for(std::map<ID, std::map <ID, std::map<ID,  
                                                std::map<ID, std::map<ID, int > > > > >::iterator iDictANs4
                                                = iDictANs3->second.begin(); 
                                                iDictANs4 != iDictANs3->second.end(); iDictANs4++)
                                            {
                                                for(std::map<ID, std::map <ID, 
                                                    std::map<ID, std::map<ID, int > > > >::iterator iDictANs5
                                                    = iDictANs4->second.begin(); 
                                                    iDictANs5 != iDictANs4->second.end(); iDictANs5++)
                                                {
                                                    for (std::map<ID, std::map<ID, std::map<ID, int > > >::iterator iDictANs6
                                                         =iDictANs5->second.begin();
                                                         iDictANs6 !=iDictANs5->second.end(); iDictANs6++)
                                                    {
                                                        for (std::map<ID, std::map<ID, int > >::iterator iDictANs7
                                                             =iDictANs6->second.begin();
                                                             iDictANs7 !=iDictANs6->second.end(); iDictANs7++)
                                                        {
                                                            for (std::map<ID, int >::iterator iDictANs8
                                                                 =iDictANs7->second.begin();
                                                                 iDictANs8 !=iDictANs7->second.end(); iDictANs8++)
                                                            {
                                                                tDictANs.push_back(allDictAngles[iDictANs8->second]);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                dLev = 2;
                                setupTargetAngleUsingdist2(tDictANs, iAN, s2, s3, dLev);
                            }
                        }
                        iAN->levelCodValue =3;
                    }
                }
                else // iFind1
                {
                    if (allDictAnglesIdx3[ha1][ha2][ha3][0].numCodValues >=10 
                                    && allDictAnglesIdx3[ha1][ha2][ha3][0].sigValue < 5.00)
                    {
                        iAN->value        = allDictAnglesIdx3[ha1][ha2][ha3][0].value;
                        iAN->sigValue     = allDictAnglesIdx3[ha1][ha2][ha3][0].sigValue;
                        iAN->numCodValues = allDictAnglesIdx3[ha1][ha2][ha3][0].numCodValues;
                    }
                    else
                    {
                        bool aCCP4S=getCCP4Angle(iAN);
                            
                        if (!aCCP4S)
                        {
                            // Finding the bond with shortest atomic distances 
                            // (or the lowest substitute costs)
                            std::vector<AngleDict> tDictANs;
                            for (std::map<ID, std::map<ID, std::map <ID, std::map<ID,  
                                 std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                                 std::map<ID, int> > > > > > > > >::iterator iDictANs1 
                                 = allDictAnglesIdx[ha1][ha2][ha3].begin();
                                 iDictANs1 !=
                                 allDictAnglesIdx[ha1][ha2][ha3].end();
                                 iDictANs1++)
                            {
                                for(std::map<ID, std::map <ID, std::map<ID, 
                                    std::map<ID, std::map<ID, std::map<ID, std::map<ID,  
                                    std::map<ID, int > > > > > > > >::iterator iDictANs2
                                    = iDictANs1->second.begin(); 
                                    iDictANs2 != iDictANs1->second.end(); iDictANs2++)
                                {
                                    for(std::map<ID, std::map <ID, std::map<ID, 
                                        std::map<ID, std::map<ID, std::map<ID,  
                                        std::map<ID, int > > > > > > >::iterator iDictANs3
                                        = iDictANs2->second.begin(); 
                                        iDictANs3 != iDictANs2->second.end(); iDictANs3++)
                                    {
                                        for(std::map<ID, std::map <ID, std::map<ID, std::map<ID, 
                                            std::map<ID, std::map<ID, int > > > > > >::iterator iDictANs4
                                            = iDictANs3->second.begin(); 
                                            iDictANs4 != iDictANs3->second.end(); iDictANs4++)
                                        {
                                            for(std::map<ID, std::map <ID, std::map<ID,  
                                                std::map<ID, std::map<ID, int > > > > >::iterator iDictANs5
                                                = iDictANs4->second.begin(); 
                                                iDictANs5 != iDictANs4->second.end(); iDictANs5++)
                                            {
                                                for(std::map<ID, std::map <ID, 
                                                    std::map<ID, std::map<ID, int > > > >::iterator iDictANs6
                                                    = iDictANs5->second.begin(); 
                                                    iDictANs6 != iDictANs5->second.end(); iDictANs6++)
                                                {
                                                    for (std::map<ID, std::map<ID, 
                                                         std::map<ID, int > > >::iterator iDictANs7
                                                         =iDictANs6->second.begin();
                                                         iDictANs7 !=iDictANs6->second.end(); iDictANs7++)
                                                    {
                                                        for (std::map<ID, std::map<ID, int > >::iterator iDictANs8
                                                             =iDictANs7->second.begin();
                                                             iDictANs8 !=iDictANs7->second.end(); iDictANs8++)
                                                        {
                                                            for (std::map<ID, int >::iterator iDictANs9
                                                                 =iDictANs8->second.begin();
                                                                 iDictANs9 !=iDictANs8->second.end(); iDictANs9++)
                                                            {
                                                                tDictANs.push_back(allDictAngles[iDictANs9->second]);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            dLev = 1;
                            setupTargetAngleUsingdist2(tDictANs, iAN, s2, s3, dLev);
                        }
                    }
                    iAN->levelCodValue =3;
                }
                
                //std::cout << "Final angle is " << iAN->value << std::endl;
                //if (iAN->hasCodValue)
                //{
                //    std::cout << "It has exact matches to a COD angle " << std::endl;
                //}
                //else
                //{
                //    std::cout << "It is generated by either averaging or shortest distance " << std::endl;
                //}
                
            }
            else   
            {
                // could not find three exact matches on 3 atomic hashing values
                // using approximate default values
                if (allAtoms[iAN->atoms[0]].bondingIdx <4)
                {
                    // std::cout << "Center atom bond index is  " << allAtoms[iAN->atoms[0]].bondingIdx<< std::endl;
                    iAN->value = DefaultOrgAngles[allAtoms[iAN->atoms[0]].bondingIdx];
                    iAN->levelCodValue = 4;
                }
                else
                {
                    std::cout << "Could not find COD angle value and even default angle value."
                        <<std::endl << "The atoms in the target angle are: " << std::endl
                        << "The inner atom "   << allAtoms[iAN->atoms[0]].id << std::endl
                        << "Its COD class is " << allAtoms[iAN->atoms[0]].codClass << std::endl
                        << "The outer atom 1 " << allAtoms[iAN->atoms[1]].id << std::endl
                        << "Its COD class is " << allAtoms[iAN->atoms[1]].codClass << std::endl
                        << "The outer atom 2 "  << allAtoms[iAN->atoms[2]].id << std::endl
                        << "Its COD class is " << allAtoms[iAN->atoms[2]].codClass << std::endl;
                
                    exit(1);
                }
                 
            }
    }
    
    void CodClassify::searchCodOrgAngles22(std::vector<AngleDict>::iterator iAN)
    {    
            /*
            codClassToAtomAng(allAtoms[iAN->atoms[0]].codClass, allAtoms[iAN->atoms[0]]);
            codClassToAtomAng(allAtoms[iAN->atoms[1]].codClass, allAtoms[iAN->atoms[1]]);
            codClassToAtomAng(allAtoms[iAN->atoms[2]].codClass, allAtoms[iAN->atoms[2]]);
            */
        
            int ha1, ha2, ha3; // s2, s3;
            ID a1NB2, a1NB, a1M, a1C, a2NB2, a2NB, a2M, a2C,a3NB2, a3NB, a3M, a3C;
            ID id2, id3;
            
            ha1   = allAtoms[iAN->atoms[0]].hashingValue;
            a1NB2 = allAtoms[iAN->atoms[0]].codNB2Symb;
            a1NB  = allAtoms[iAN->atoms[0]].codNBSymb;
            a1M   = allAtoms[iAN->atoms[0]].codAtmMain;
            a1C   = allAtoms[iAN->atoms[0]].codClass;
                     
            if ((int)allAtoms[iAN->atoms[1]].hashingValue
                    <(int)allAtoms[iAN->atoms[2]].hashingValue )
            {
                //s2 =1;
                //s3 =2;
                
                ha2  =allAtoms[iAN->atoms[1]].hashingValue;
                ha3  =allAtoms[iAN->atoms[2]].hashingValue;
           
                a2NB2= allAtoms[iAN->atoms[1]].codNB2Symb;
                a3NB2= allAtoms[iAN->atoms[2]].codNB2Symb;
          
                a2NB = allAtoms[iAN->atoms[1]].codNBSymb;
                a3NB = allAtoms[iAN->atoms[2]].codNBSymb;
         
                a2M  = allAtoms[iAN->atoms[1]].codAtmMain;
                a3M  = allAtoms[iAN->atoms[2]].codAtmMain;
                
                a2C  = allAtoms[iAN->atoms[1]].codClass;
                a3C  = allAtoms[iAN->atoms[2]].codClass;
                
                id2  = allAtoms[iAN->atoms[1]].id;
                id3  = allAtoms[iAN->atoms[2]].id;
            }
            else if ((int)allAtoms[iAN->atoms[1]].hashingValue
                    ==(int)allAtoms[iAN->atoms[2]].hashingValue)
            {
                if (allAtoms[iAN->atoms[1]].codAtmMain.size() 
                     > allAtoms[iAN->atoms[2]].codAtmMain.size())
                {
                    ha2  =allAtoms[iAN->atoms[1]].hashingValue;
                    ha3  =allAtoms[iAN->atoms[2]].hashingValue;
           
                    a2NB2= allAtoms[iAN->atoms[1]].codNB2Symb;
                    a3NB2= allAtoms[iAN->atoms[2]].codNB2Symb;
          
                    a2NB = allAtoms[iAN->atoms[1]].codNBSymb;
                    a3NB = allAtoms[iAN->atoms[2]].codNBSymb;
         
                    a2M  = allAtoms[iAN->atoms[1]].codAtmMain;
                    a3M  = allAtoms[iAN->atoms[2]].codAtmMain;
                
                    a2C  = allAtoms[iAN->atoms[1]].codClass;
                    a3C  = allAtoms[iAN->atoms[2]].codClass;
                
                    id2  = allAtoms[iAN->atoms[1]].id;
                    id3  = allAtoms[iAN->atoms[2]].id;
                    
                }
                else if (allAtoms[iAN->atoms[1]].codAtmMain.size() 
                     <  allAtoms[iAN->atoms[2]].codAtmMain.size())
                {
                    ha2  =allAtoms[iAN->atoms[2]].hashingValue;
                    ha3  =allAtoms[iAN->atoms[1]].hashingValue;
           
                    a2NB2= allAtoms[iAN->atoms[2]].codNB2Symb;
                    a3NB2= allAtoms[iAN->atoms[1]].codNB2Symb;
          
                    a2NB = allAtoms[iAN->atoms[2]].codNBSymb;
                    a3NB = allAtoms[iAN->atoms[1]].codNBSymb;
         
                    a2C  = allAtoms[iAN->atoms[2]].codClass;
                    a3C  = allAtoms[iAN->atoms[1]].codClass;
                
                    id2  = allAtoms[iAN->atoms[2]].id;
                    id3  = allAtoms[iAN->atoms[1]].id;
                }
                else
                {
                    std::vector<sortMap4> vM;
                
                    sortMap4   m1, m2;
                
                    m1.id    = allAtoms[iAN->atoms[1]].id;
                    m2.id    = allAtoms[iAN->atoms[2]].id;
                    m1.ha    =(int)allAtoms[iAN->atoms[1]].hashingValue;
                    m2.ha    =(int)allAtoms[iAN->atoms[2]].hashingValue;
                    m1.lev2  = allAtoms[iAN->atoms[1]].codNB2Symb;
                    m2.lev2  = allAtoms[iAN->atoms[2]].codNB2Symb;
                    m1.lev3  = allAtoms[iAN->atoms[1]].codNBSymb;
                    m2.lev3  = allAtoms[iAN->atoms[2]].codNBSymb;
                    m1.key   = allAtoms[iAN->atoms[1]].codAtmMain;
                    m2.key   = allAtoms[iAN->atoms[2]].codAtmMain;
                    m1.lev4  = allAtoms[iAN->atoms[1]].codClass;
                    m2.lev4  = allAtoms[iAN->atoms[2]].codClass;
                
                    vM.push_back(m1);
                    vM.push_back(m2);
                
                    std::sort(vM.begin(), vM.end(), sortMapkey4);
                
                    id2    = vM[0].id;
                    id3    = vM[1].id;
                    ha2    = vM[0].ha;
                    ha3    = vM[1].ha;
                    a2NB2  = vM[0].lev2;
                    a3NB2  = vM[1].lev2;
                    a2NB   = vM[0].lev3;
                    a3NB   = vM[1].lev3;
                    a2M    = vM[0].key;
                    a3M    = vM[1].key;
                    a2C    = vM[0].lev4;
                    a3C    = vM[1].lev4;
                }
            }
            else
            {
                //s2 =2;
                //s3 =1;
                ha2  =allAtoms[iAN->atoms[2]].hashingValue;
                ha3  =allAtoms[iAN->atoms[1]].hashingValue;
           
                a2NB2= allAtoms[iAN->atoms[2]].codNB2Symb;
                a3NB2= allAtoms[iAN->atoms[1]].codNB2Symb;
          
                a2NB = allAtoms[iAN->atoms[2]].codNBSymb;
                a3NB = allAtoms[iAN->atoms[1]].codNBSymb;
         
                a2M  = allAtoms[iAN->atoms[2]].codAtmMain;
                a3M  = allAtoms[iAN->atoms[1]].codAtmMain;
                
                a2C  = allAtoms[iAN->atoms[2]].codClass;
                a3C  = allAtoms[iAN->atoms[1]].codClass;
                
                id2  = allAtoms[iAN->atoms[2]].id;
                id3  = allAtoms[iAN->atoms[1]].id;
                
            }
            
            std::cout << "\n%=========================================%" << std::endl;
            std::cout << "Angle between " << allAtoms[iAN->atoms[0]].id 
                      << "(center) and " << id2 << " and " 
                      << id3 << std::endl << std::endl;
            
            std::cout << "Center Atom : " <<  allAtoms[iAN->atoms[0]].id << std::endl
                      <<  "Its Hashing code  " << ha1 <<std::endl
                      <<  "Its codNBSymb " <<  a1NB << " its codNB2Symb " << a1NB2 << std::endl
                      <<  "atom type main sec " << a1M << std::endl
                      <<  "atom type full " << a1C << std::endl<< std::endl;
            
            std::cout << "Atom2 : "  << id2 << std::endl
                      << "Its Hashing code "  << ha2    << std::endl
                      << "Its codNBSymb " <<  a2NB  << " its codNB2Symb " << a2NB2 << std::endl
                      << "Its atom type main section " << a2M << std::endl
                      << " Its atom type full " <<  a2C << std::endl << std::endl; 
            
            std::cout << "Atom3: "   << id3 << std::endl
                      << "Its Hashing code  " <<  ha3   <<std::endl
                      << "Its codNBSymb " <<  a3NB  << " its codNB2Symb "  << a3NB2 << std::endl
                      << "Its atom type main section " << a3M << std::endl
                      << "Its atom type full " << a3C << std::endl << std::endl;
            
            // int dLev = 0;
            
                    
            if(allDictAnglesIdxD.find(ha1) != allDictAnglesIdxD.end() &&
               allDictAnglesIdxD[ha1].find(ha2) != allDictAnglesIdxD[ha1].end() &&
               allDictAnglesIdxD[ha1][ha2].find(ha3) != allDictAnglesIdxD[ha1][ha2].end())
            {
                    
                //std::map<ID,  std::map<ID,  std::map<ID, 
                //std::vector<aValueSet> > > >::iterator  iFind1 
                //=allDictAnglesIdx2[ha1][ha2][ha3].find(a1NB2);
                
                
                if(allDictAnglesIdxD[ha1][ha2][ha3].find(a1NB2) != allDictAnglesIdxD[ha1][ha2][ha3].end())
                {
                    // iFind1 a1NB2 matches
                    std::cout << "Found " << a1NB2 << std::endl;
                    // a1NB2 
                  
                    if ( allDictAnglesIdxD[ha1][ha2][ha3][a1NB2].find(a2NB2)!=allDictAnglesIdxD[ha1][ha2][ha3][a1NB2].end())
                    {
                        //iFind2: a1NB2, a2NB2 match
                        std::cout << "Found " << a2NB2 << std::endl;
                        
                        if (allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2].find(a3NB2)
                            != allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2].end())
                        {
                            //iFind3: a1NB2, a2NB2 and a3NB2 match
                            std::cout << "Found " << a3NB2 << std::endl;
                            if ( allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2].find(a1NB) 
                                 !=allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2].end())
                            {
     
                                // iFind4: a1NB2, a2NB2, a3NB2, a1NB match
                                std::cout << "Found " << a1NB << std::endl;
                                
                                if ( allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB].find(a2NB)
                                     != allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB].end())
                                {
                                    // iFind5: a1NB2, a2NB2, a3NB2, a1NB, a2NB match
                                    std::cout << "Found " << a2NB << std::endl;
                                   
                                    if ( allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB].find(a3NB)
                                         !=allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB].end())
                                    {
                                        //iFind6 a1NB2, a2NB2, a3NB2, a1NB, a2NB, a3NB match
                                        std::cout << "Found " << a3NB << std::endl;
                                        
                                        if(allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB].find(a1M)
                                            != allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB].end())
                                        {
                                            // iFind7: a1NB2, a2NB2, a3NB2, a1NB, a2NB, a3NB, a1M match
                                            std::cout << "Found " << a1M << std::endl;
                                            if(allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1M].find(a2M)
                                               != allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1M].end())
                                            {
                                                //iFind8: a1NB2, a2NB2, a3NB2, a1NB, a2NB, a3NB, a1M, a2M match
                                                std::cout << "Found " << a2M << std::endl;
                                                
                                                if( allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1M][a2M].find(a3M)!=
                                                    allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1M][a2M].end())
                                                {
                                                    //iFind9: a1NB2, a2NB2, a3NB2, a1NB, a2NB, a3NB, a1M, a2M and a3M all match
                                                    std::cout << "Found " << a3M << std::endl;
                                                    if( allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1M][a2M][a3M].find(a1C)!=
                                                        allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1M][a2M][a3M].end())
                                                    {
                                                        //iFind10: a1NB2, a2NB2, a3NB2, a1NB, a2NB, a3NB, a1M, a2M, a3M, a1C all match
                                                        std::cout << "Found " << a1C << std::endl;
                                                        
                                                        if(allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1M][a2M][a3M][a1C].find(a2C)!=
                                                           allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1M][a2M][a3M][a1C].end())
                                                        {
                                                            //iFind11: a1NB2, a2NB2, a3NB2, a1NB, a2NB, a3NB, a1M, a2M, a3M, a1C, a2C all match
                                                            std::cout << "Found " << a2C << std::endl;
                                                            if(allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1M][a2M][a3M][a1C][a2C].find(a3C)!=
                                                               allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][a1M][a2M][a3M][a1C][a2C].end())
                                                            {
                                                                //iFind12:  all keys match
                                                                std::cout << "Found " << a3C << std::endl;
                                                                // COD has such an angle value
                                                                int iPos = allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2]
                                                                           [a1NB][a2NB][a3NB][a1M][a2M][a3M][a1C][a2C][a3C];
                                                                //std::cout << "iPos =" << iPos << std::endl;
                                                                if (allDictAnglesD[iPos].numCodValues > 5)
                                                                {
                                                                    iAN->value        = allDictAnglesD[iPos].value;
                                                                    iAN->sigValue     = allDictAnglesD[iPos].sigValue;
                                                                    iAN->numCodValues = allDictAnglesD[iPos].numCodValues;
                                                                    iAN->hasCodValue  = true;
                                                                    iAN->levelCodValue = 0;
                                                                    std::cout << "COD finds the exact value " << iAN->value << std::endl;
                                                                    std::cout << "the sigma value " << iAN->sigValue << std::endl;
                                                                }
                                                                else
                                                                {
                                                                    iAN->value = allDictAnglesIdx1D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2]
                                                                                 [a1NB][a2NB][a3NB][a1M][a2M][a3M][0].value;
                                                                    iAN->sigValue = allDictAnglesIdx1D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2]
                                                                                    [a1NB][a2NB][a3NB][a1M][a2M][a3M][0].sigValue;
                                                                    iAN->numCodValues = allDictAnglesIdx1D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2]
                                                                                        [a1NB][a2NB][a3NB][a1M][a2M][a3M][0].numCodValues;
                                                                    iAN->hasCodValue  = true;
                                                                    iAN->levelCodValue = 1;
                                                                }
                                                            }
                                                            else
                                                            {
                                                                // iFind12 failed a3C
                                                                iAN->value = allDictAnglesIdx1D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2]
                                                                                 [a1NB][a2NB][a3NB][a1M][a2M][a3M][0].value;
                                                                iAN->sigValue = allDictAnglesIdx1D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2]
                                                                                 [a1NB][a2NB][a3NB][a1M][a2M][a3M][0].sigValue;
                                                                iAN->numCodValues = allDictAnglesIdx1D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2]
                                                                                 [a1NB][a2NB][a3NB][a1M][a2M][a3M][0].numCodValues;
                                                                iAN->hasCodValue  = true;
                                                                iAN->levelCodValue = 1;   
                                                            }
                                                        }
                                                        else
                                                        {
                                                            //iFind11 failed a2C
                                                            iAN->value = allDictAnglesIdx1D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2]
                                                                         [a1NB][a2NB][a3NB][a1M][a2M][a3M][0].value;
                                                            iAN->sigValue = allDictAnglesIdx1D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2]
                                                                         [a1NB][a2NB][a3NB][a1M][a2M][a3M][0].sigValue;
                                                            iAN->numCodValues = allDictAnglesIdx1D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2]
                                                                         [a1NB][a2NB][a3NB][a1M][a2M][a3M][0].numCodValues;
                                                            iAN->hasCodValue  = true;
                                                            iAN->levelCodValue = 1;
                                                        }                  
                                                    }
                                                    else
                                                    {
                                                        //iFind10 failed a1C
                                                        iAN->value = allDictAnglesIdx1D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2]
                                                                                       [a1NB][a2NB][a3NB][a1M][a2M][a3M][0].value;
                                                        iAN->sigValue = allDictAnglesIdx1D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2]
                                                                                          [a1NB][a2NB][a3NB][a1M][a2M][a3M][0].sigValue;
                                                        iAN->numCodValues = allDictAnglesIdx1D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2]
                                                                                        [a1NB][a2NB][a3NB][a1M][a2M][a3M][0].numCodValues;
                                                        iAN->hasCodValue  = true;
                                                        iAN->levelCodValue = 1;
                                                        
                                                    }    
                                                }
                                                else // iFind9 failed a3M
                                                {
                                                    //aValueSet   tVaS;
                                                    //setValueSet(tVaS, allDictAnglesIdx2D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB]);
                                                    
                                                    iAN->value         = allDictAnglesIdx2D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][0].value;
                                                    iAN->sigValue      = allDictAnglesIdx2D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][0].sigValue;
                                                    if (iAN->sigValue < 0.1 || iAN->sigValue > 3.0)
                                                    {
                                                       iAN->sigValue = 3.0;
                                                    }
                                                    iAN->numCodValues  = allDictAnglesIdx2D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][0].numCodValues;
                                                    iAN->levelCodValue = 1;
                                                    
                                                }
                                            }
                                            else // iFind8 failed a2M
                                            {
                                                //aValueSet   tVaS;
                                                //setValueSet(tVaS, allDictAnglesIdx2D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB]);
                                                    
                                                iAN->value        = allDictAnglesIdx2D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][0].value;        
                                                iAN->sigValue     = allDictAnglesIdx2D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][0].sigValue;
                                                if (iAN->sigValue < 0.1 || iAN->sigValue > 3.0)
                                                {
                                                    iAN->sigValue = 3.0;
                                                }
                                                iAN->numCodValues = allDictAnglesIdx2D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][0].numCodValues;
                                                    
                                                iAN->levelCodValue = 1;
                                                
                                            }
                                        }
                                        else //iFind7 failed a1M
                                        {
                                            
                                            //aValueSet   tVaS;
                                            //setValueSet(tVaS, allDictAnglesIdx2D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB]);
                                                    
                                            iAN->value          = allDictAnglesIdx2D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][0].value;        
                                            iAN->sigValue       = allDictAnglesIdx2D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][0].sigValue;
                                            if (iAN->sigValue < 0.1 || iAN->sigValue > 3.0)
                                            {
                                                iAN->sigValue = 3.0;
                                            }
                                            iAN->numCodValues   = allDictAnglesIdx2D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB][a3NB][0].numCodValues;                        
                                            iAN->levelCodValue = 1;
                                            
                                        }
                                        
                                    }
                                    else // iFind6 failed a3NB 
                                    {
                                        /*
                                        std::vector<aValueSet> tDictANs;
                                        for (std::map<ID, std::vector<aValueSet> >::iterator iDictANs1 
                                             = allDictAnglesIdx2D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB].begin();
                                             iDictANs1 !=
                                             allDictAnglesIdx2D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB][a2NB].end();
                                             iDictANs1++)
                                        {
                                            for(std::vector<aValueSet>::iterator iDictANs2
                                                = iDictANs1->second.begin(); 
                                                iDictANs2 != iDictANs1->second.end(); iDictANs2++)
                                            {
                                                tDictANs.push_back(*iDictANs2);
                                            }
                                        }

                                        aValueSet   tVaS;
                                            
                                        setValueSet(tVaS, tDictANs);
                                        **/
                                        
                                        iAN->value             = allDictAnglesIdx3D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][0].value;        
                                        iAN->sigValue          = allDictAnglesIdx3D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][0].sigValue;
                                        if (iAN->sigValue < 0.1 || iAN->sigValue > 3.0)
                                        {
                                            iAN->sigValue = 3.0;
                                        }
                                        iAN->numCodValues      = allDictAnglesIdx3D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][0].numCodValues;                        
                                        
                                        iAN->levelCodValue = 2;
                                    }
                                }
                                else // iFind5 failed a2NB
                                {   /*
                                    std::vector<aValueSet> tDictANs;
                                    
                                    for (std::map<ID, std::map<ID, std::vector<aValueSet> > >::iterator iDictANs1 
                                             = allDictAnglesIdx2D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB].begin();
                                             iDictANs1 !=
                                             allDictAnglesIdx2D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][a1NB].end();
                                             iDictANs1++)
                                    {
                                        for(std::map<ID, std::vector<aValueSet> >::iterator iDictANs2
                                                = iDictANs1->second.begin(); 
                                                iDictANs2 != iDictANs1->second.end(); iDictANs2++)
                                        {
                                            for(std::vector<aValueSet>::iterator iDictANs3
                                                   = iDictANs2->second.begin(); 
                                                  iDictANs3 != iDictANs2->second.end(); iDictANs3++)
                                            {
                                                tDictANs.push_back(*iDictANs2);
                                            }
                                            
                                        }
                                    }

                                    aValueSet   tVaS;
                                            
                                    setValueSet(tVaS, tDictANs);
                                    */                
                                    iAN->value             = allDictAnglesIdx3D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][0].value;        
                                    iAN->sigValue          = allDictAnglesIdx3D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][0].sigValue;
                                    if (iAN->sigValue < 0.1 || iAN->sigValue > 3.0)
                                    {
                                        iAN->sigValue = 3.0;
                                    }
                                    iAN->numCodValues      = allDictAnglesIdx3D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][0].numCodValues;                           
                                    iAN->levelCodValue = 2;
                                }
                            }
                            else // iFind4 a1NB
                            {
                                /*
                                std::vector<aValueSet> tDictANs;
                                    
                                for (std::map<ID, std::map<ID, std::map<ID, std::vector<aValueSet> > > >::iterator iDictANs1 
                                       = allDictAnglesIdx2D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2].begin();
                                     iDictANs1 != allDictAnglesIdx2D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2].end();
                                             iDictANs1++)
                                {
                                    for(std::map<ID, std::map<ID, std::vector<aValueSet> > >::iterator iDictANs2
                                            = iDictANs1->second.begin(); 
                                            iDictANs2 != iDictANs1->second.end(); iDictANs2++)
                                    {
                                        for(std::map<ID, std::vector<aValueSet> >::iterator iDictANs3
                                             = iDictANs2->second.begin(); 
                                            iDictANs3 != iDictANs2->second.end(); iDictANs3++)
                                        { 
                                            for(std::vector<aValueSet>::iterator iDictANs4 =
                                                  iDictANs3->second.begin();
                                                iDictANs4 != iDictANs3->second.end(); iDictANs4++)
                                            {
                                                tDictANs.push_back(*iDictANs4);
                                            }
                                        }   
                                    }
                                }

                                aValueSet   tVaS;
                                            
                                setValueSet(tVaS, tDictANs);
                                */
                                
                                iAN->value             = allDictAnglesIdx3D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][0].value;        
                                iAN->sigValue          = allDictAnglesIdx3D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][0].sigValue;
                                if (iAN->sigValue < 0.1 || iAN->sigValue > 3.0)
                                {
                                    iAN->sigValue  = 3.0;
                                }
                                iAN->numCodValues  = allDictAnglesIdx3D[ha1][ha2][ha3][a1NB2][a2NB2][a3NB2][0].numCodValues;                           
                                    
                                iAN->levelCodValue = 2;
                                
                            }
                        }
                        else // iFind3 failed a3NB2
                        {
                            
                            std::vector<aValueSet> tDictANs;
                                    
                            for (std::map<ID, std::map<ID, std::map<ID, std::map<ID,
                                 std::map<ID, std::map<ID, std::map<ID, std::map<ID,
                                 std::map<ID, std::map<ID, int> > > > > > > > > >::iterator iDictANs1 
                                       = allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2].begin();
                                   iDictANs1 != allDictAnglesIdxD[ha1][ha2][ha3][a1NB2][a2NB2].end();
                                             iDictANs1++)
                            {
                                for (std::map<ID, std::map<ID, std::map<ID, std::map<ID,
                                     std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                                     std::map<ID, int> > > > > > > > >::iterator iDictANs2 
                                     = iDictANs1->second.begin(); iDictANs2 !=iDictANs1->second.end(); iDictANs2++)
                                {
                                    for (std::map<ID, std::map<ID, std::map<ID,
                                         std::map<ID, std::map<ID, std::map<ID,
                                         std::map<ID, std::map<ID, int> > > > > > > >::iterator iDictANs3 
                                         = iDictANs2->second.begin(); iDictANs3 !=iDictANs2->second.end(); iDictANs3++)
                                    {
                                        for (std::map<ID, std::map<ID, std::map<ID,
                                             std::map<ID, std::map<ID, std::map<ID,
                                             std::map<ID, int> > > > > > >::iterator iDictANs4 
                                             = iDictANs3->second.begin(); iDictANs4 !=iDictANs3->second.end(); iDictANs4++)
                                        {
                                            for (std::map<ID, std::map<ID, std::map<ID,
                                                 std::map<ID, std::map<ID, 
                                                 std::map<ID, int> > > > > >::iterator iDictANs5 
                                                 = iDictANs4->second.begin(); iDictANs5 !=iDictANs4->second.end(); iDictANs5++)
                                            {
                                                for (std::map<ID, std::map<ID, std::map<ID, std::map<ID,  
                                                     std::map<ID, int> > > > >::iterator iDictANs6 
                                                 = iDictANs5->second.begin(); iDictANs6 !=iDictANs5->second.end(); iDictANs6++)
                                                {
                                                    for (std::map<ID, std::map<ID, std::map<ID,   
                                                         std::map<ID, int> > > >::iterator iDictANs7 
                                                         = iDictANs6->second.begin(); 
                                                         iDictANs7 !=iDictANs6->second.end(); iDictANs7++)
                                                    {
                                                       for (std::map<ID, std::map<ID, std::map<ID, int> > >::iterator iDictANs8 
                                                            = iDictANs7->second.begin(); 
                                                            iDictANs8 !=iDictANs7->second.end(); iDictANs8++)
                                                       {
                                                           for (std::map<ID, std::map<ID, int> >::iterator iDictANs9 
                                                            = iDictANs8->second.begin(); 
                                                            iDictANs9 !=iDictANs8->second.end(); iDictANs9++)
                                                           {
                                                               for (std::map<ID, int>::iterator iDictANs10 
                                                               = iDictANs9->second.begin(); 
                                                               iDictANs10 !=iDictANs9->second.end(); iDictANs10++)
                                                               {
                                                                   tDictANs.push_back(allDictAnglesD[iDictANs10->second]);
                                                               }
                                                           }
                                                       }      
                                                    }
                                                }
                                            }
                                            
                                         }    
                                    }
                                    
                                }
                            }
                                    

                            aValueSet   tVaS;
                                            
                            setValueSet(tVaS, tDictANs);
                            
                            if (tVaS.numCodValues >=10 && tVaS.sigValue < 5.00)
                            {
                                iAN->value        = tVaS.value;
                                iAN->sigValue     = tVaS.sigValue;
                                iAN->numCodValues = tVaS.numCodValues;
                            }
                            else
                            {
                                bool aCCP4S=getCCP4Angle(iAN);
                            
                                if (!aCCP4S)
                                {
                                    iAN->value        = tVaS.value;
                                    iAN->sigValue     = tVaS.sigValue;
                                    iAN->numCodValues = tVaS.numCodValues;
                                }
                            }
                            iAN->levelCodValue = 3;
                        }
                    }
                    else // iFind2 a2NB2
                    {
                        std::cout << "iFind2 " <<  a2NB2 << std::endl;
                        std::vector<aValueSet> tDictANs;
                                    
                        for (std::map<ID, std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                             std::map<ID, std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                             std::map<ID, int> > > > > > > > > > >::iterator iDictANs1 
                                = allDictAnglesIdxD[ha1][ha2][ha3][a1NB2].begin();
                             iDictANs1 != allDictAnglesIdxD[ha1][ha2][ha3][a1NB2].end();
                             iDictANs1++)
                        {
                            for (std::map<ID, std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                                 std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                                 std::map<ID, int> > > > > > > > > >::iterator iDictANs2 
                                 = iDictANs1->second.begin(); iDictANs2 !=iDictANs1->second.end(); iDictANs2++)
                                {
                                    for (std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                                         std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                                         std::map<ID, int> > > > > > > > >::iterator iDictANs3 
                                         = iDictANs2->second.begin(); iDictANs3 !=iDictANs2->second.end(); iDictANs3++)
                                    {
                                        for (std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                                             std::map<ID, std::map<ID, std::map<ID,
                                             std::map<ID, int> > > > > > > >::iterator iDictANs4 
                                             = iDictANs3->second.begin(); iDictANs4 !=iDictANs3->second.end(); iDictANs4++)
                                        {
                                            for (std::map<ID, std::map<ID, std::map<ID,
                                                 std::map<ID, std::map<ID, std::map<ID, 
                                                 std::map<ID, int> > > > > > >::iterator iDictANs5 
                                                 = iDictANs4->second.begin(); iDictANs5 !=iDictANs4->second.end(); iDictANs5++)
                                            {
                                                for (std::map<ID, std::map<ID, std::map<ID, std::map<ID,  
                                                     std::map<ID, std::map<ID, int> > > > > >::iterator iDictANs6 
                                                 = iDictANs5->second.begin(); iDictANs6 !=iDictANs5->second.end(); iDictANs6++)
                                                {
                                                    for (std::map<ID, std::map<ID, std::map<ID, std::map<ID,   
                                                         std::map<ID, int> > > > >::iterator iDictANs7 
                                                         = iDictANs6->second.begin(); 
                                                         iDictANs7 !=iDictANs6->second.end(); iDictANs7++)
                                                    {
                                                       for (std::map<ID, std::map<ID, 
                                                            std::map<ID, std::map<ID, int> > > >::iterator iDictANs8 
                                                            = iDictANs7->second.begin(); 
                                                            iDictANs8 !=iDictANs7->second.end(); iDictANs8++)
                                                       {
                                                           for (std::map<ID, std::map<ID, std::map<ID, int> > >::iterator iDictANs9 
                                                            = iDictANs8->second.begin(); 
                                                            iDictANs9 !=iDictANs8->second.end(); iDictANs9++)
                                                           {
                                                               for (std::map<ID, std::map<ID, int> >::iterator iDictANs10 
                                                               = iDictANs9->second.begin(); 
                                                               iDictANs10 !=iDictANs9->second.end(); iDictANs10++)
                                                               {
                                                                   for (std::map<ID, int>::iterator iDictANs11 
                                                                        = iDictANs10->second.begin(); 
                                                                        iDictANs11 !=iDictANs10->second.end(); iDictANs11++)
                                                                   {
                                                                        tDictANs.push_back(allDictAnglesD[iDictANs11->second]);
                                                                   }
                                                               }
                                                           }
                                                       }      
                                                    }
                                                }
                                            }
                                            
                                         }    
                                    }
                                    
                                }
                            }
                                    

                            aValueSet   tVaS;
                            setValueSet(tVaS, tDictANs);
                            
                            if (tVaS.numCodValues >=10 
                                 && tVaS.sigValue < 5.00)
                            {
                                iAN->value        = tVaS.value;
                                iAN->sigValue     = tVaS.sigValue;
                                iAN->numCodValues = tVaS.numCodValues;
                            }
                            else
                            {
                                bool aCCP4S=getCCP4Angle(iAN);
                            
                                if (!aCCP4S)
                                {
                                    iAN->value        = tVaS.value;
                                    iAN->sigValue     = tVaS.sigValue;
                                    iAN->numCodValues = tVaS.numCodValues;
                                    
                                    // Finding the bond with shortest atomic distances 
                                    // (or the lowest substitute costs)
                                    /*
                                    std::vector<AngleDict> tDictANs;
                                    for (std::map<ID, std::map<ID, std::map <ID, std::map<ID,  
                                     std::map<ID, std::map<ID, std::map<ID, 
                                     std::map<ID, int> > > > > > > >::iterator iDictANs1 
                                     = allDictAnglesIdx[ha1][ha2][ha3][a1NB2].begin();
                                     iDictANs1 !=
                                     allDictAnglesIdx[ha1][ha2][ha3][a1NB2].end();
                                     iDictANs1++)
                                {
                                    for(std::map<ID, std::map <ID, std::map<ID, 
                                        std::map<ID, std::map<ID, std::map<ID,  
                                        std::map<ID, int > > > > > > >::iterator iDictANs2
                                        = iDictANs1->second.begin(); 
                                        iDictANs2 != iDictANs1->second.end(); iDictANs2++)
                                    {
                                        for(std::map<ID, std::map <ID, std::map<ID, std::map<ID, 
                                            std::map<ID, std::map<ID, int > > > > > >::iterator iDictANs3
                                            = iDictANs2->second.begin(); 
                                            iDictANs3 != iDictANs2->second.end(); iDictANs3++)
                                        {
                                            for(std::map<ID, std::map <ID, std::map<ID,  
                                                std::map<ID, std::map<ID, int > > > > >::iterator iDictANs4
                                                = iDictANs3->second.begin(); 
                                                iDictANs4 != iDictANs3->second.end(); iDictANs4++)
                                            {
                                                for(std::map<ID, std::map <ID, 
                                                    std::map<ID, std::map<ID, int > > > >::iterator iDictANs5
                                                    = iDictANs4->second.begin(); 
                                                    iDictANs5 != iDictANs4->second.end(); iDictANs5++)
                                                {
                                                    for (std::map<ID, std::map<ID, std::map<ID, int > > >::iterator iDictANs6
                                                         =iDictANs5->second.begin();
                                                         iDictANs6 !=iDictANs5->second.end(); iDictANs6++)
                                                    {
                                                        for (std::map<ID, std::map<ID, int > >::iterator iDictANs7
                                                             =iDictANs6->second.begin();
                                                             iDictANs7 !=iDictANs6->second.end(); iDictANs7++)
                                                        {
                                                            for (std::map<ID, int >::iterator iDictANs8
                                                                 =iDictANs7->second.begin();
                                                                 iDictANs8 !=iDictANs7->second.end(); iDictANs8++)
                                                            {
                                                                tDictANs.push_back(allDictAngles[iDictANs8->second]);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                                dLev = 2;
                                setupTargetAngleUsingdist2(tDictANs, iAN, s2, s3, dLev);
                                */
                            }
                        }
                        iAN->levelCodValue =3;
                    }
                }
                else // iFind1 a1NB2
                {
                    //if (allDictAnglesIdx3[ha1][ha2][ha3][0].numCodValues >=10 
                    //                && allDictAnglesIdx3[ha1][ha2][ha3][0].sigValue < 5.00)
                   // {
                   //     iAN->value        = allDictAnglesIdx3[ha1][ha2][ha3][0].value;
                   //     iAN->sigValue     = allDictAnglesIdx3[ha1][ha2][ha3][0].sigValue;
                   //     iAN->numCodValues = allDictAnglesIdx3[ha1][ha2][ha3][0].numCodValues;
                   // }
                   // else
                   //
                    bool aCCP4S=getCCP4Angle(iAN);
                            
                    if (!aCCP4S)
                    {
                        std::vector<aValueSet> tDictANs;
                                    
                        for (std::map<ID, std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                             std::map<ID, std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                             std::map<ID, std::map<ID, int> > > > > > > > > > > >::iterator iDictANs1 
                                = allDictAnglesIdxD[ha1][ha2][ha3].begin();
                             iDictANs1 != allDictAnglesIdxD[ha1][ha2][ha3].end();
                             iDictANs1++)
                        {
                            for (std::map<ID, std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                                 std::map<ID, std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                                 std::map<ID, int> > > > > > > > > > >::iterator iDictANs2 
                                 = iDictANs1->second.begin(); iDictANs2 !=iDictANs1->second.end(); iDictANs2++)
                                {
                                    for (std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                                         std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                                         std::map<ID, std::map<ID, int> > > > > > > > > >::iterator iDictANs3 
                                         = iDictANs2->second.begin(); iDictANs3 !=iDictANs2->second.end(); iDictANs3++)
                                    {
                                        for (std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                                             std::map<ID, std::map<ID, std::map<ID, std::map<ID,
                                             std::map<ID, int> > > > > > > > >::iterator iDictANs4 
                                             = iDictANs3->second.begin(); iDictANs4 !=iDictANs3->second.end(); iDictANs4++)
                                        {
                                            for (std::map<ID, std::map<ID, std::map<ID,
                                                 std::map<ID, std::map<ID, std::map<ID, std::map<ID,  
                                                 std::map<ID, int> > > > > > > >::iterator iDictANs5 
                                                 = iDictANs4->second.begin(); iDictANs5 !=iDictANs4->second.end(); iDictANs5++)
                                            {
                                                for (std::map<ID, std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                                                     std::map<ID, std::map<ID, int> > > > > > >::iterator iDictANs6 
                                                 = iDictANs5->second.begin(); iDictANs6 !=iDictANs5->second.end(); iDictANs6++)
                                                {
                                                    for (std::map<ID, std::map<ID, std::map<ID, std::map<ID,   
                                                         std::map<ID,std::map<ID, int> > > > > >::iterator iDictANs7 
                                                         = iDictANs6->second.begin(); 
                                                         iDictANs7 !=iDictANs6->second.end(); iDictANs7++)
                                                    {
                                                       for (std::map<ID, std::map<ID, std::map<ID,  
                                                            std::map<ID, std::map<ID, int> > > > >::iterator iDictANs8 
                                                            = iDictANs7->second.begin(); 
                                                            iDictANs8 !=iDictANs7->second.end(); iDictANs8++)
                                                       {
                                                           for (std::map<ID, std::map<ID, std::map<ID,  
                                                                std::map<ID, int> > > >::iterator iDictANs9 
                                                            = iDictANs8->second.begin(); 
                                                            iDictANs9 !=iDictANs8->second.end(); iDictANs9++)
                                                           {
                                                               for (std::map<ID, std::map<ID,  
                                                                    std::map<ID, int> > >::iterator iDictANs10 
                                                               = iDictANs9->second.begin(); 
                                                               iDictANs10 !=iDictANs9->second.end(); iDictANs10++)
                                                               {
                                                                   for (std::map<ID, std::map<ID, int> >::iterator iDictANs11 
                                                                        = iDictANs10->second.begin(); 
                                                                        iDictANs11 !=iDictANs10->second.end(); iDictANs11++)
                                                                   {
                                                                       for (std::map<ID, int>::iterator iDictANs12 
                                                                            = iDictANs11->second.begin(); 
                                                                            iDictANs12 !=iDictANs11->second.end(); iDictANs12++)
                                                                       {
                                                                           tDictANs.push_back(allDictAnglesD[iDictANs12->second]);
                                                                       }
                                                                   }
                                                               }
                                                           }
                                                       }      
                                                    }
                                                }
                                            }
                                            
                                         }    
                                    }
                                    
                                }
                            }
                                    

                            aValueSet   tVaS;
                                            
                            setValueSet(tVaS, tDictANs);
                            
                            iAN->value        = tVaS.value;
                            iAN->sigValue     = tVaS.sigValue;
                            iAN->numCodValues = tVaS.numCodValues;
                            
                            // Finding the bond with shortest atomic distances 
                            // (or the lowest substitute costs)
                            
                            /*
                            std::vector<AngleDict> tDictANs;
                            for (std::map<ID, std::map<ID, std::map <ID, std::map<ID,  
                                 std::map<ID, std::map<ID, std::map<ID, std::map<ID, 
                                 std::map<ID, int> > > > > > > > >::iterator iDictANs1 
                                 = allDictAnglesIdx[ha1][ha2][ha3].begin();
                                 iDictANs1 !=
                                 allDictAnglesIdx[ha1][ha2][ha3].end();
                                 iDictANs1++)
                            {
                                for(std::map<ID, std::map <ID, std::map<ID, 
                                    std::map<ID, std::map<ID, std::map<ID, std::map<ID,  
                                    std::map<ID, int > > > > > > > >::iterator iDictANs2
                                    = iDictANs1->second.begin(); 
                                    iDictANs2 != iDictANs1->second.end(); iDictANs2++)
                                {
                                    for(std::map<ID, std::map <ID, std::map<ID, 
                                        std::map<ID, std::map<ID, std::map<ID,  
                                        std::map<ID, int > > > > > > >::iterator iDictANs3
                                        = iDictANs2->second.begin(); 
                                        iDictANs3 != iDictANs2->second.end(); iDictANs3++)
                                    {
                                        for(std::map<ID, std::map <ID, std::map<ID, std::map<ID, 
                                            std::map<ID, std::map<ID, int > > > > > >::iterator iDictANs4
                                            = iDictANs3->second.begin(); 
                                            iDictANs4 != iDictANs3->second.end(); iDictANs4++)
                                        {
                                            for(std::map<ID, std::map <ID, std::map<ID,  
                                                std::map<ID, std::map<ID, int > > > > >::iterator iDictANs5
                                                = iDictANs4->second.begin(); 
                                                iDictANs5 != iDictANs4->second.end(); iDictANs5++)
                                            {
                                                for(std::map<ID, std::map <ID, 
                                                    std::map<ID, std::map<ID, int > > > >::iterator iDictANs6
                                                    = iDictANs5->second.begin(); 
                                                    iDictANs6 != iDictANs5->second.end(); iDictANs6++)
                                                {
                                                    for (std::map<ID, std::map<ID, 
                                                         std::map<ID, int > > >::iterator iDictANs7
                                                         =iDictANs6->second.begin();
                                                         iDictANs7 !=iDictANs6->second.end(); iDictANs7++)
                                                    {
                                                        for (std::map<ID, std::map<ID, int > >::iterator iDictANs8
                                                             =iDictANs7->second.begin();
                                                             iDictANs8 !=iDictANs7->second.end(); iDictANs8++)
                                                        {
                                                            for (std::map<ID, int >::iterator iDictANs9
                                                                 =iDictANs8->second.begin();
                                                                 iDictANs9 !=iDictANs8->second.end(); iDictANs9++)
                                                            {
                                                                tDictANs.push_back(allDictAngles[iDictANs9->second]);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                              
                            }
                            dLev = 1;
                            setupTargetAngleUsingdist2(tDictANs, iAN, s2, s3, dLev);
                             */
                        }
                    
                    iAN->levelCodValue =3;
                
                }
                //std::cout << "Final angle is " << iAN->value << std::endl;
                //if (iAN->hasCodValue)
                //{
                //    std::cout << "It has exact matches to a COD angle " << std::endl;
                //}
                //else
                //{
                //    std::cout << "It is generated by either averaging or shortest distance " << std::endl;
                //}
                
            }
            else   
            {
                // could not find three exact matches on 3 atomic hashing values
                // using approximate default values
                if (allAtoms[iAN->atoms[0]].bondingIdx <4)
                {
                    // std::cout << "Center atom bond index is  " << allAtoms[iAN->atoms[0]].bondingIdx<< std::endl;
                    iAN->value = DefaultOrgAngles[allAtoms[iAN->atoms[0]].bondingIdx];
                    iAN->levelCodValue = 4;
                }
                else
                {
                    std::cout << "Could not find COD angle value and even default angle value."
                        <<std::endl << "The atoms in the target angle are: " << std::endl
                        << "The inner atom "   << allAtoms[iAN->atoms[0]].id << std::endl
                        << "Its COD class is " << allAtoms[iAN->atoms[0]].codClass << std::endl
                        << "The outer atom 1 " << allAtoms[iAN->atoms[1]].id << std::endl
                        << "Its COD class is " << allAtoms[iAN->atoms[1]].codClass << std::endl
                        << "The outer atom 2 "  << allAtoms[iAN->atoms[2]].id << std::endl
                        << "Its COD class is " << allAtoms[iAN->atoms[2]].codClass << std::endl;
                
                    exit(1);
                }
                 
            }
    }
    
    
    
    bool CodClassify::getCCP4Angle(std::vector<AngleDict>::iterator iAN)
    {
        bool aFound = false;
        
        ID atmCen = allAtoms[iAN->atoms[0]].id;
        std::cout << "center atom CCP4 type  " << atmCen << std::endl; 
        if (ccp4Angles.find(atmCen) !=ccp4Angles.end())
        {
            ID atm1 = allAtoms[iAN->atoms[1]].id;
            ID atm2 = allAtoms[iAN->atoms[2]].id;
             
            if (ccp4Angles[atmCen].find(atm1) !=ccp4Angles[atmCen].end())
            {
                std::cout << "atom 1 " << atm1;
                if (ccp4Angles[atmCen][atm1].find(atm2) 
                        != ccp4Angles[atmCen][atm1].end() )
                {
                    std::cout << " atom 2 " << atm2 << std::endl;
                    iAN->value = ccp4Angles[atmCen][atm1][atm2];
                    aFound = true;
                }
            }
            else if (ccp4Angles[atmCen].find(atm2) !=ccp4Angles[atmCen].end())
            {
                std::cout << "atom 1 " << atm2;
                if (ccp4Angles[atmCen][atm2].find(atm1) 
                        != ccp4Angles[atmCen][atm2].end() )
                {
                    std::cout << " atom 2 " << atm1 << std::endl;
                    iAN->value = ccp4Angles[atmCen][atm2][atm1];
                    aFound = true;
                }
            }      
        }
        
        return aFound;
    }
    
    /*
    void CodClassify::searchOneOrgAngleFromCodUsingSqlite(sqlite3* tCombDB, 
                                                          std::vector<AngleDict>::iterator iAN)
    {
            int ha1, ha2, ha3;
            ID a1NB2, a1NB, a1C, a2NB2, a2NB, a2C,a3NB2, a3NB, a3C;
            ID id2, id3;
            ha1   = allAtoms[iAN->atoms[0]].hashingValue;
            a1NB2 = allAtoms[iAN->atoms[0]].codNB2Symb;
            a1NB  = allAtoms[iAN->atoms[0]].codNBSymb;
            a1C   = allAtoms[iAN->atoms[0]].codClass;
            
            std::cout << "Atom1 " <<  allAtoms[iAN->atoms[0]].id 
                      << " Its Cod class " <<  a1C <<  " and Hashing "
                      << ha1 <<std::endl;
            std::cout << " Its codNBSymb " <<  a1NB << " its codNB2Symb " << a1NB2 << std::endl;
           
            
            if ((int)allAtoms[iAN->atoms[1]].hashingValue
                    <(int)allAtoms[iAN->atoms[2]].hashingValue )
            {
                ha2  =allAtoms[iAN->atoms[1]].hashingValue;
                ha3  =allAtoms[iAN->atoms[2]].hashingValue;
           
                a2NB2= allAtoms[iAN->atoms[1]].codNB2Symb;
                a3NB2= allAtoms[iAN->atoms[2]].codNB2Symb;
          
                a2NB = allAtoms[iAN->atoms[1]].codNBSymb;
                a3NB = allAtoms[iAN->atoms[2]].codNBSymb;
         
                a2C  = allAtoms[iAN->atoms[1]].codClass;
                a3C  = allAtoms[iAN->atoms[2]].codClass;
                
                id2  = allAtoms[iAN->atoms[1]].id;
                id3  = allAtoms[iAN->atoms[2]].id;
            }
            else if ((int)allAtoms[iAN->atoms[1]].hashingValue
                    ==(int)allAtoms[iAN->atoms[2]].hashingValue)
            {
                if ((int)allAtoms[iAN->atoms[1]].codClass.size() <=
                        (int)allAtoms[iAN->atoms[2]].codClass.size())
                {
                    ha2  =allAtoms[iAN->atoms[1]].hashingValue;
                    ha3  =allAtoms[iAN->atoms[2]].hashingValue;
           
                    a2NB2= allAtoms[iAN->atoms[1]].codNB2Symb;
                    a3NB2= allAtoms[iAN->atoms[2]].codNB2Symb;
          
                    a2NB = allAtoms[iAN->atoms[1]].codNBSymb;
                    a3NB = allAtoms[iAN->atoms[2]].codNBSymb;
         
                    a2C  = allAtoms[iAN->atoms[1]].codClass;
                    a3C  = allAtoms[iAN->atoms[2]].codClass;
                
                    id2  = allAtoms[iAN->atoms[1]].id;
                    id3  = allAtoms[iAN->atoms[2]].id;
                }
                else
                {
                    ha2  =allAtoms[iAN->atoms[2]].hashingValue;
                    ha3  =allAtoms[iAN->atoms[1]].hashingValue;
           
                    a2NB2= allAtoms[iAN->atoms[2]].codNB2Symb;
                    a3NB2= allAtoms[iAN->atoms[1]].codNB2Symb;
          
                    a2NB = allAtoms[iAN->atoms[2]].codNBSymb;
                    a3NB = allAtoms[iAN->atoms[1]].codNBSymb;
         
                    a2C  = allAtoms[iAN->atoms[2]].codClass;
                    a3C  = allAtoms[iAN->atoms[1]].codClass;
                
                    id2  = allAtoms[iAN->atoms[2]].id;
                    id3  = allAtoms[iAN->atoms[1]].id;
                }
            }
            else
            {
                ha2  =allAtoms[iAN->atoms[2]].hashingValue;
                ha3  =allAtoms[iAN->atoms[1]].hashingValue;
           
                a2NB2= allAtoms[iAN->atoms[2]].codNB2Symb;
                a3NB2= allAtoms[iAN->atoms[1]].codNB2Symb;
          
                a2NB = allAtoms[iAN->atoms[2]].codNBSymb;
                a3NB = allAtoms[iAN->atoms[1]].codNBSymb;
         
                a2C  = allAtoms[iAN->atoms[2]].codClass;
                a3C  = allAtoms[iAN->atoms[1]].codClass;
                
                id2  = allAtoms[iAN->atoms[2]].id;
                id3  = allAtoms[iAN->atoms[1]].id;
                
            }
            
            
            std::cout << "Atom2 " <<  id2 
                      << " Its Cod class " <<  a2C << " and Hashing " << ha2 <<std::endl;
            std::cout << " Its codNBSymb " <<  a2NB << " its codNB2Symb " << a2NB2 << std::endl;
            
            std::cout << "Atom3 " <<  id3 
                      << " Its Cod class " <<  a3C  << " and Hashing "    << ha3 <<std::endl;
            std::cout << " Its codNBSymb " <<  a3NB << " its codNB2Symb " << a3NB2 << std::endl;
            
            
            int dLev = 0;
            
            std::map<ID, ID>  propNB;
            propNB["a1NB2"] = a1NB2;
            propNB["a2NB2"] = a2NB2;
            propNB["a3NB2"] = a3NB2;
            propNB["a1NB"]  = a1NB;
            propNB["a2NB"]  = a2NB;
            propNB["a3NB"]  = a3NB;
            propNB["ha1"]   = IntToStr(ha1);
            propNB["ha2"]   = IntToStr(ha2);
            propNB["ha3"]   = IntToStr(ha3);
            
            // Now query the databases
        
            // 1. search for exact match of atom classes.
            std::string qQue = "SELECT * from angles WHERE atomCen=\'" 
                              + a1C + "\' and atom1=\'" + a2C
                              + "\' and atom2=\'" + a3C
                              + "\';"; // one query statement 
            // std::cout << "1st qQue is " << qQue << std::endl;
        
        
            std::vector<std::vector<std::string> > qResults;
        
            sqlite3Query(tCombDB, qQue.c_str(), qResults);
            
            if ((int)qResults.size() !=0)
            {
                // find exact match
                int nVals= StrToInt(qResults[0][14]);
                if (nVals >=4)
                {
                    iAN->valueST      = StrToReal(qResults[0][12]);
                    iAN->value        = iAN->valueST;
                    iAN->sigValueST   = StrToReal(qResults[0][13]);
                    iAN->sigValue     = iAN->sigValueST;
                    iAN->numCodValues = nVals;
                }
                else
                {
                    dLev =8;
                    searchOneOrgAngleUsingSqliteL(tCombDB, iAN, dLev, propNB, a1C); 
                }
            }
            else
            {
                dLev = 8;
                searchOneOrgAngleUsingSqliteL(tCombDB, iAN, dLev, propNB, a1C);
            }               
    }
    
    bool CodClassify::searchOneOrgAngleUsingSqliteL(sqlite3* tCombDB, 
                                                    std::vector<AngleDict>::iterator tA, 
                                                    int& tLev, std::map<ID,ID> tPropNB,
                                                    std::string  tA1C)
    {
        bool iFind = false;
        std::string tQue;
        std::vector<std::vector<std::string> > qResults;
        
        // set different query strings for different levels of searching
        if (tLev >0)
        {
            setQueStrAn(tQue, tLev, tPropNB, tA1C);
            sqlite3Query(tCombDB, tQue.c_str(), qResults);
            if ((int)qResults.size() !=0)
            { 
                setOneAngleByMean(qResults, tA);
                iFind = true;
            }
            else
            {
               tLev -=1;
               iFind = searchOneOrgAngleUsingSqliteL(tCombDB,  tA, tLev,
                                                     tPropNB, tA1C);
            }
        }
          
        return iFind;
        
    }
    */
    
    void CodClassify::setQueStrAn(std::string& tQue, int tLev, 
                                  std::map<ID,ID> tPropNB, std::string tA1C)
    {
        if (tLev==8)
        {
            tQue = "SELECT * from angles WHERE hash1=" 
                    + tPropNB["ha1"] + " and hash2=" + tPropNB["ha2"] + " and hash3=" + tPropNB["ha3"]
                    + " and NeighB1N=\'" + tPropNB["a1NB2"] + "\' and NeighB2N=\'" + tPropNB["a2NB2"]
                    + "\' and NeighB3N=\'" + tPropNB["a3NB2"] + "\' and NeighB1C=\'" + tPropNB["a1NB"] 
                    +"\' and NeighB2C=\'"  + tPropNB["a2NB"]  + "\' and NeighB3C=\'" + tPropNB["a3NB"] 
                    + "\' and atomCen=\'"    + tA1C + "\';";
        }
        else if (tLev==7)
        {
            tQue = "SELECT * from angles WHERE hash1=" 
                    + tPropNB["ha1"] + " and hash2=" + tPropNB["ha2"] + " and hash3=" + tPropNB["ha3"]
                    + " and NeighB1N=\'" + tPropNB["a1NB2"] + "\' and NeighB2N=\'" + tPropNB["a2NB2"]
                    + "\' and NeighB3N=\'" + tPropNB["a3NB2"] + "\' and NeighB1C=\'" + tPropNB["a1NB"] 
                    +"\' and NeighB2C=\'" + tPropNB["a2NB"] + "\' and NeighB3C=\'" + tPropNB["a3NB"] 
                    + "\';";
        }
        else if (tLev==6)
        {
            tQue = "SELECT * from angles WHERE hash1=" 
                    + tPropNB["ha1"] + " and hash2=" + tPropNB["ha2"] + " and hash3=" + tPropNB["ha3"]
                    + " and NeighB1N=\'" + tPropNB["a1NB2"] + "\' and NeighB2N=\'" + tPropNB["a2NB2"]
                    + "\' and NeighB3N=\'" + tPropNB["a3NB2"] + "\' and NeighB1C=\'" + tPropNB["a1NB"] 
                    +"\' and NeighB2C=\'" + tPropNB["a2NB"] + "\';";
        }
        else if(tLev==5)
        {
            tQue = "SELECT * from angles WHERE hash1=" 
                    + tPropNB["ha1"] + " and hash2=" + tPropNB["ha2"] + " and hash3=" + tPropNB["ha3"]
                    + " and NeighB1N=\'"   + tPropNB["a1NB2"] + "\' and NeighB2N=\'" + tPropNB["a2NB2"]
                    + "\' and NeighB3N=\'" + tPropNB["a3NB2"] + "\' and NeighB1C=\'" + tPropNB["a1NB"] 
                    + "\';";
        }
        else if (tLev==4)
        {
            tQue = "SELECT * from angles WHERE hash1=" 
                    + tPropNB["ha1"] + " and hash2=" + tPropNB["ha2"] + " and hash3=" + tPropNB["ha3"]
                    + " and NeighB1N=\'" + tPropNB["a1NB2"] + "\' and NeighB2N=\'" + tPropNB["a2NB2"]
                    + "\' and NeighB3N=\'" + tPropNB["a3NB2"]+ "\';";
        }
        else if(tLev==3)
        {
            tQue = "SELECT * from angles WHERE hash1=" 
                    + tPropNB["ha1"] + " and hash2=" + tPropNB["ha2"] + " and hash3=" + tPropNB["ha3"]
                    + " and NeighB1N=\'" + tPropNB["a1NB2"] + "\' and NeighB2N=\'" + tPropNB["a2NB2"]
                    + "\';";
        }
        else if(tLev==2)
        {
            tQue = "SELECT * from angles WHERE hash1=" 
                    + tPropNB["ha1"] + " and hash2=" + tPropNB["ha2"] + " and hash3=" + tPropNB["ha3"]
                    + " and NeighB1N=\'" + tPropNB["a1NB2"] 
                    + "\';";
        }
        else if(tLev==1)
        {
            tQue = "SELECT * from angles WHERE hash1=" 
                    + tPropNB["ha1"] + " and hash2=" + tPropNB["ha2"] + " and hash3=" + tPropNB["ha3"]
                    + ";";
        }
        
    }
    
    void CodClassify::setOneAngleByMean(std::vector<std::vector<std::string> >& tResults, 
                                        std::vector<AngleDict>::iterator tA)
    {
        int nAngles=0;
        REAL sumVals=0.0;
        for (int i=0; i < (int)tResults.size(); i++)
        {
            int n1       = StrToInt(tResults[i][14]); 
            nAngles  += n1;
            sumVals +=(StrToReal(tResults[i][12])*n1);
        }
        
        if (nAngles >0)
        {
            REAL tMean = sumVals/nAngles;
            tA->valueST = tMean;
            tA->value   = tMean;
            REAL tCur  = 0.0;
            for (int i=0; i < (int)tResults.size(); i++)
            {
                int n1  = StrToInt(tResults[i][14]);
                REAL m1 = StrToReal(tResults[i][12]);
                for (int j=0; j < n1; j++)
                {
                    tCur +=((tMean-m1)*(tMean-m1));
                }
            }
            
            if (nAngles ==1)
            {
                tA->sigValue = sqrt(tCur);
            }
            else
            {
                tA->sigValueST = sqrt(tCur/(nAngles-1));
            }
            tA->sigValue = tA->sigValueST;
        }

    }
    
    void CodClassify::searchCodAnglesWithNonCenteredMetal(std::vector<AngleDict>::iterator iAN)
    {
        ID id1, idNB1;
        id1= allAtoms[iAN->atoms[0]].chemType;
        
        ID tIdNB1, tIdNB2, tIdNB3;
        // the following is temporal scheme. write a function for it if the scheme adopted
        int tM =0;
        for (std::vector<int>::iterator iNB =allAtoms[iAN->atoms[0]].connAtoms.begin();
                iNB != allAtoms[iAN->atoms[0]].connAtoms.end(); iNB++)
        {
            if (allAtoms[*iNB].isMetal)
            {
                tM++;
            }
        }
        idNB1 = "m"+IntToStr(tM) + "n" + IntToStr((int)allAtoms[iAN->atoms[0]].connAtoms.size()-tM);

        tM =0;
        for (std::vector<int>::iterator iNB =allAtoms[iAN->atoms[1]].connAtoms.begin();
                iNB != allAtoms[iAN->atoms[1]].connAtoms.end(); iNB++)
        {
            if (allAtoms[*iNB].isMetal)
            {
                tM++;
            }
        }
        tIdNB2 = "m" + IntToStr(tM) + "n" + IntToStr((int)allAtoms[iAN->atoms[1]].connAtoms.size()-tM);
        
        tM =0;
        for (std::vector<int>::iterator iNB =allAtoms[iAN->atoms[2]].connAtoms.begin();
                iNB != allAtoms[iAN->atoms[2]].connAtoms.end(); iNB++)
        {
            if (allAtoms[*iNB].isMetal)
            {
                tM++;
            }
        }
        tIdNB3 = "m" + IntToStr(tM) + "n" + IntToStr((int)allAtoms[iAN->atoms[2]].connAtoms.size()-tM);        
        
        ID id2, idNB2, id3, idNB3;
        std::list<ID> tIds;
        tIds.push_back(allAtoms[iAN->atoms[1]].chemType);
        tIds.push_back(allAtoms[iAN->atoms[2]].chemType);
        //std::cout << "initial front " << tIds.front() << std::endl;
        //std::cout << "back " << tIds.back() << std::endl;
        tIds.sort();
        //std::cout << "After sorting: front " << tIds.front() << std::endl;
        //std::cout << "back " << tIds.back() << std::endl;
        if(tIds.front() == allAtoms[iAN->atoms[1]].chemType)
        {
            id2   = allAtoms[iAN->atoms[1]].chemType;
            idNB2 = tIdNB2; 
            id3   = allAtoms[iAN->atoms[2]].chemType;
            idNB3 = tIdNB3;
        }
        else
        {
            id2 = allAtoms[iAN->atoms[2]].chemType;
            idNB2 = tIdNB3;
            id3 = allAtoms[iAN->atoms[1]].chemType;
            idNB3 = tIdNB2;
        }
        
        std::cout << "Center atom ID " << id1 << " NB: " << idNB1 <<std::endl;
        std::cout << "outer atom 1 ID " << id2 << " NB: " << idNB2 <<std::endl;
        std::cout << "outer atom 2 ID " << id3 << " NB: " << idNB3 <<std::endl;
                    
        std::map<ID,  std::map<ID,  std::map<ID, 
        std::map<ID,  std::map<ID,  std::map<ID,
        int > > > > > >::iterator  iFind1 
                 =allDictNonCenMetAnglesIdx.find(id1);
        if (iFind1 !=allDictNonCenMetAnglesIdx.end())
        {
            std::map<ID,  std::map<ID,  std::map<ID, 
            std::map<ID,  std::map<ID, int > > > > >::iterator  iFind2 
                 =allDictNonCenMetAnglesIdx[id1].find(idNB1);
            if (iFind2 !=allDictNonCenMetAnglesIdx[id1].end())
            {
                std::map<ID,  std::map<ID,  std::map<ID, 
                std::map<ID,  int > > > >::iterator  iFind3 
                 =allDictNonCenMetAnglesIdx[id1][idNB1].find(id2);
                if (iFind3 !=allDictNonCenMetAnglesIdx[id1][idNB1].end())
                {
                    std::map<ID,  std::map<ID,  std::map<ID, 
                    int > > >::iterator  iFind4 
                     =allDictNonCenMetAnglesIdx[id1][idNB1][id2].find(idNB2);
                    if(iFind4 !=allDictNonCenMetAnglesIdx[id1][idNB1][id2].end())
                    {
                        std::map<ID,  std::map<ID,  int > >::iterator  iFind5 
                          =allDictNonCenMetAnglesIdx[id1][idNB1][id2][idNB2].find(id3);
                        if (iFind5 !=allDictNonCenMetAnglesIdx[id1][idNB1][id2][idNB2].end())
                        {
                            std::map<ID, int>::iterator  iFind6 
                            =allDictNonCenMetAnglesIdx[id1][idNB1][id2][idNB2][id3].find(idNB3);
                            if (iFind6 !=allDictNonCenMetAnglesIdx[id1][idNB1][id2][idNB2][id3].end())
                            {
                                int ida = allDictNonCenMetAnglesIdx[id1][idNB1][id2][idNB2][id3][idNB3];
                                if (allDictAngles[ida].numCodValues > 5)
                                {
                                    iAN->value        = allDictAngles[ida].value;
                                    iAN->sigValue     = allDictAngles[ida].sigValue;
                                    iAN->hasCodValue  = true;
                                    iAN->numCodValues = allDictAngles[ida].numCodValues;
                                }
                                else
                                { 
                                    std::vector<ID>  tAngsKeys;
                                    std::vector<int> tAngsIdx;
                                    for (std::map<ID, int>::iterator iA=allDictNonCenMetAnglesIdx[id1][idNB1][id2][idNB2][id3].begin();
                                            iA != allDictNonCenMetAnglesIdx[id1][idNB1][id2][idNB2][id3].end(); iA++)
                                    {
                                        tAngsKeys.push_back(iA->first);
                                        tAngsIdx.push_back(iA->second);
                                    }
                                    ID tKey = getMatchedKey(tAngsKeys, idNB3, tAngsIdx);
                                    if ((int)tKey.size() ==4)
                                    {
                                        int ida    = allDictNonCenMetAnglesIdx[id1][idNB1][id2][idNB2][id3][tKey];
                                        iAN->value    = allDictAngles[ida].value;
                                        iAN->sigValue = allDictAngles[ida].sigValue;
                                        iAN->hasCodValue  = true;
                                        iAN->numCodValues = allDictAngles[ida].numCodValues;
                                    }
                                    else
                                    {
                                        // Find no similar keys
                                        std::vector<AngleDict> tAngs;
                                        for (std::vector<int>::iterator iN=tAngsIdx.begin();
                                                iN !=tAngsIdx.end(); iN++)
                                        {
                                            tAngs.push_back(allDictAngles[*iN]);
                                        }
                                        setupTargetAngleUsingMean(tAngs, iAN);
                                    }
                                }
                            
                            }
                            else
                            {
                                // No match for idNB3
                                
                                std::vector<ID>  tAngsKeys;
                                std::vector<int> tAngsIdx;
                                for (std::map<ID, int>::iterator iA=allDictNonCenMetAnglesIdx[id1][idNB1][id2][idNB2][id3].begin();
                                            iA != allDictNonCenMetAnglesIdx[id1][idNB1][id2][idNB2][id3].end(); iA++)
                                {
                                    tAngsKeys.push_back(iA->first);
                                    tAngsIdx.push_back(iA->second);
                                }
                                ID tKey = getMatchedKey(tAngsKeys, idNB3, tAngsIdx);
                                if ((int)tKey.size() ==4)
                                {
                                    int ida    = allDictNonCenMetAnglesIdx[id1][idNB1][id2][idNB2][id3][tKey];
                                    iAN->value    = allDictAngles[ida].value;
                                    iAN->sigValue = allDictAngles[ida].sigValue;
                                    iAN->hasCodValue  = true;
                                    iAN->numCodValues = allDictAngles[ida].numCodValues;
                                }
                                else
                                {
                                    // Find no similar keys
                                    std::vector<AngleDict> tAngs;
                                    for (std::vector<int>::iterator iN=tAngsIdx.begin();
                                           iN !=tAngsIdx.end(); iN++)
                                    {
                                        tAngs.push_back(allDictAngles[*iN]);
                                    }
                                    setupTargetAngleUsingMean(tAngs, iAN);
                                }
                            }
                        }
                        else
                        {
                            // Not find id3
                            std::vector<std::vector<ID> > tKeyPairs;
                            std::vector<int> tAngIdxs;
                            std::vector<ID> targetPair;
                            targetPair.push_back(id3);
                            targetPair.push_back(idNB3);
                            int resultIdx =-1;
                            
                            for (std::map<ID,  std::map<ID,  int > >::iterator 
                                 iA=allDictNonCenMetAnglesIdx[id1][idNB1][id2][idNB2].begin();
                                 iA!=allDictNonCenMetAnglesIdx[id1][idNB1][id2][idNB2].end();
                                 iA++)
                            {
                                for (std::map<ID, int>::iterator iA1=iA->second.begin();
                                        iA1 != iA->second.end(); iA1++)
                                {
                                    std::vector<ID> aKeyPair;
                                    aKeyPair.push_back(iA->first);
                                    aKeyPair.push_back(iA1->first);
                                    tKeyPairs.push_back(aKeyPair);
                                    tAngIdxs.push_back(iA1->second);
                                }
                            }
                            
                            getMatched2dVect(tKeyPairs, targetPair, tAngIdxs, resultIdx); 
                            
                            if (resultIdx !=-1)
                            {
                                iAN->value =allDictAngles[resultIdx].value;
                                iAN->value =allDictAngles[resultIdx].sigValue;
                                iAN->hasCodValue = false;
                                iAN->numCodValues = allDictAngles[resultIdx].numCodValues;
                            }
                            else
                            {
                                // using the mean value
                                std::vector<AngleDict> tAngs;
                                for (std::vector<int>::iterator iN=tAngIdxs.begin();
                                           iN !=tAngIdxs.end(); iN++)
                                {
                                    tAngs.push_back(allDictAngles[*iN]);
                                }
                                setupTargetAngleUsingMean(tAngs, iAN);
                            }
                        }
                    }
                    else
                    {
                        // No match for idNB2
                        std::vector<std::vector<ID> > tKeyPairs;
                        std::vector<int> tAngIdxs;
                        std::vector<ID> targetPair;
                        targetPair.push_back(idNB2);
                        targetPair.push_back(id3);
                        targetPair.push_back(idNB3);
                        int resultIdx =-1;
                        for (std::map<ID, std::map<ID,  std::map<ID,  int > > >::iterator
                                iA = allDictNonCenMetAnglesIdx[id1][idNB1][id2].begin();
                                iA != allDictNonCenMetAnglesIdx[id1][idNB1][id2].end();
                                iA++)
                        {
                            for (std::map<ID,  std::map<ID,  int > >::iterator 
                                 iA1 =iA->second.begin(); iA1!=iA->second.end(); iA1++)
                            {
                                for (std::map<ID, int>::iterator iA2=iA1->second.begin();
                                        iA2 != iA1->second.end(); iA2++)
                                {
                                    std::vector<ID> aKeyPair;
                                    aKeyPair.push_back(iA->first);
                                    aKeyPair.push_back(iA1->first);
                                    aKeyPair.push_back(iA2->first);
                                    tKeyPairs.push_back(aKeyPair);
                                    tAngIdxs.push_back(iA2->second);
                                }
                            }
                        }
                        
                        getMatched3dVect(tKeyPairs, targetPair, tAngIdxs, resultIdx);
                        if (resultIdx !=-1)
                        {
                            iAN->value =allDictAngles[resultIdx].value;
                            iAN->value =allDictAngles[resultIdx].sigValue;
                            iAN->hasCodValue = false;
                            iAN->numCodValues = allDictAngles[resultIdx].numCodValues;
                        }
                        else
                        {
                            // using the mean value
                            std::vector<AngleDict> tAngs;
                            for (std::vector<int>::iterator iN=tAngIdxs.begin();
                                           iN !=tAngIdxs.end(); iN++)
                            {
                                tAngs.push_back(allDictAngles[*iN]);
                            }
                            setupTargetAngleUsingMean(tAngs, iAN);
                        }
                    }
                }
                else
                {
                    // not find id2
                    std::vector<std::vector<ID> > tKeyPairs;
                    std::vector<int> tAngIdxs;
                    std::vector<ID> targetPair;
                    targetPair.push_back(id2);
                    targetPair.push_back(idNB2);
                    targetPair.push_back(id3);
                    targetPair.push_back(idNB3);
                    int resultIdx =-1;
                    
                    for (std::map<ID, std::map<ID, std::map<ID,  std::map<ID,  int > > > >::iterator
                         iA = allDictNonCenMetAnglesIdx[id1][idNB1].begin();
                         iA != allDictNonCenMetAnglesIdx[id1][idNB1].end(); iA++)
                    {
                        for (std::map<ID, std::map<ID,  std::map<ID,  int > > >::iterator 
                                 iA1 =iA->second.begin(); iA1!=iA->second.end(); iA1++)
                        {
                            for (std::map<ID,  std::map<ID,  int > >::iterator 
                                 iA2 =iA1->second.begin(); iA2!=iA1->second.end();
                                 iA2++)
                            {
                                for (std::map<ID, int>::iterator iA3=iA2->second.begin();
                                        iA3 != iA2->second.end(); iA3++)
                                {
                                    std::vector<ID> aKeyPair;
                                    aKeyPair.push_back(iA->first);
                                    aKeyPair.push_back(iA1->first);
                                    aKeyPair.push_back(iA2->first);
                                    aKeyPair.push_back(iA3->first);
                                    tKeyPairs.push_back(aKeyPair);
                                    tAngIdxs.push_back(iA3->second);
                                }
                            }
                        }
                    }
                        
                    getMatched4dVect(tKeyPairs, targetPair, tAngIdxs, resultIdx);
                    if (resultIdx !=-1)
                    {
                        iAN->value =allDictAngles[resultIdx].value;
                        iAN->value =allDictAngles[resultIdx].sigValue;
                        iAN->hasCodValue = false;
                        iAN->numCodValues = allDictAngles[resultIdx].numCodValues;
                    }
                    else
                    {
                        // using the mean value
                        std::vector<AngleDict> tAngs;
                        for (std::vector<int>::iterator iN=tAngIdxs.begin();
                             iN !=tAngIdxs.end(); iN++)
                        {
                            tAngs.push_back(allDictAngles[*iN]);
                        }
                        setupTargetAngleUsingMean(tAngs, iAN);
                    }
                }
            }
            else
            {
                // could not find idNB1
                std::vector<std::vector<ID> > tKeyPairs;
                std::vector<int> tAngIdxs;
                std::vector<ID> targetPair;
                targetPair.push_back(idNB1);
                targetPair.push_back(id2);
                targetPair.push_back(idNB2);
                targetPair.push_back(id3);
                targetPair.push_back(idNB3);
                int resultIdx =-1;
                    
                for (std::map<ID, std::map<ID, std::map<ID,  
                     std::map<ID, std::map<ID,  int > > > > >::iterator
                     iA = allDictNonCenMetAnglesIdx[id1].begin();
                     iA != allDictNonCenMetAnglesIdx[id1].end(); iA++)
                {
                    for (std::map<ID, std::map<ID,  std::map<ID,  
                         std::map<ID, int > > > >::iterator 
                         iA1 =iA->second.begin(); iA1!=iA->second.end(); iA1++)
                    {
                        for (std::map<ID,  std::map<ID, std::map<ID,  int > > >::iterator 
                             iA2 =iA1->second.begin(); iA2!=iA1->second.end(); iA2++)
                        {
                            for (std::map<ID, std::map<ID, int> >::iterator iA3=iA2->second.begin();
                                 iA3 != iA2->second.end(); iA3++)
                            {
                                for (std::map<ID, int>::iterator iA4=iA3->second.begin();
                                iA4 != iA3->second.end(); iA4++)
                                {
                                    std::vector<ID> aKeyPair;
                                    aKeyPair.push_back(iA->first);
                                    aKeyPair.push_back(iA1->first);
                                    aKeyPair.push_back(iA2->first);
                                    aKeyPair.push_back(iA3->first);
                                    aKeyPair.push_back(iA4->first);
                                    tKeyPairs.push_back(aKeyPair);
                                    tAngIdxs.push_back(iA4->second);
                                }
                            }
                        }
                    }
                }
                        
               getMatched5dVect(tKeyPairs, targetPair, tAngIdxs, resultIdx);
               if (resultIdx !=-1)
               {
                   iAN->value =allDictAngles[resultIdx].value;
                   iAN->value =allDictAngles[resultIdx].sigValue;
                   iAN->hasCodValue = false;
                   iAN->numCodValues = allDictAngles[resultIdx].numCodValues;
               }
               else
               {
                   // using the mean value
                   std::vector<AngleDict> tAngs;
                   for (std::vector<int>::iterator iN=tAngIdxs.begin();
                         iN !=tAngIdxs.end(); iN++)
                   {
                       tAngs.push_back(allDictAngles[*iN]);
                   }
                   setupTargetAngleUsingMean(tAngs, iAN);
               }
            }
                
        }
        else
        {
            std::cout << "Not finding the centered organic element " << id1 
                    << " in the dictionary of the angles with non-centered metal elements."
                    << std::endl << "That is impossible !" << std::endl;
            exit(1);
        }
    }
    
    
    // The following is used finding a key with the same number of non-metal NB,
    // the closest number of metal NB
    ID  CodClassify::getMatchedKey(std::vector<ID> tKeys, ID tTarget, std::vector<int> tIdxs)
    {
        int tM = StrToInt(tTarget.substr(1,1)); // metal neighbors
        ID  tN = tTarget.substr(3,1); // non-metal neighbors
        
        ID  fID="";    // final return ID
        int tMin = 100, tM1=-1;
        for (int i=0; i < (int)tKeys.size(); i++)
        {
            if (tKeys[i].substr(3,1)==tN && allDictAngles[tIdxs[i]].numCodValues > 5)
            {
                int tM2 = StrToInt(tKeys[i].substr(1,1));
                int tAbs = abs(tM2-tM);
                if (tAbs < tMin)
                {
                    tMin = tAbs;
                    tM1  = tM2;
                }
            }
        }
        if (tM1 !=-1)
        {
            fID = "m" + IntToStr(tM1) + "n"+tN;
        }
        
        return fID;
    }
    
    void CodClassify::getMatched2dVect(std::vector<std::vector<ID> > tPairs, std::vector<ID> aPair, 
                                       std::vector<int> tIdxs, int & rIdx)
    {
        int minScore  = 1000000;
        
        int groupUnit = 100;
        int nonMUnit  = 50;
        int metUnit   = 20;
        
        PeriodicTable tP;
        
        for (int i=0; i < (int)tPairs.size(); i++)
        {
            int tScore = 0;
            //  group difference
            tScore += (abs(tP.elements[tPairs[i][0]]["group"]-tP.elements[aPair[0]]["group"])*groupUnit);
            int n1  = StrToInt(tPairs[i][1].substr(3,1));
            int n2  = StrToInt(aPair[1].substr(3,1));
            tScore += (abs(n1-n2)*nonMUnit);
            int m1 = StrToInt(tPairs[i][1].substr(1,1));
            int m2 = StrToInt(aPair[1].substr(1,1));
            tScore += (abs(m1-m2)*metUnit);
            int tCNum = allDictAngles[tIdxs[i]].numCodValues;
            if (tScore < minScore && tCNum > 5)
            {
                minScore = tScore;
                rIdx     = i;  
            }
        }
    }
   
    
    void CodClassify::getMatched3dVect(std::vector<std::vector<ID> > tPairs, 
            std::vector<ID> aPair, std::vector<int> tIdxs, int& rIdx)
    {
            
        int minScore  = 1000000;
        
        int groupUnit = 100;
        int nonMUnit  = 50;
        int metUnit   = 20;
        
        PeriodicTable tP;
        
        for (int i=0; i < (int)tPairs.size(); i++)
        {
            int tScore = 0;
            int n1  = StrToInt(tPairs[i][0].substr(3,1));
            int n2  = StrToInt(aPair[0].substr(3,1));
            tScore += (abs(n1-n2)*nonMUnit);
            int m1  = StrToInt(tPairs[i][0].substr(1,1));
            int m2  = StrToInt(aPair[0].substr(1,1));
            tScore += (abs(m1-m2)*metUnit);
            tScore += (abs(tP.elements[tPairs[i][1]]["group"]-tP.elements[aPair[1]]["group"])*groupUnit);
            int n3  = StrToInt(tPairs[i][2].substr(3,1));
            int n4  = StrToInt(aPair[2].substr(3,1));
            tScore += (abs(n3-n4)*nonMUnit);
            int m3  = StrToInt(tPairs[i][2].substr(1,1));
            int m4  = StrToInt(aPair[2].substr(1,1));
            tScore += (abs(m3-m4)*metUnit);
            int tCNum = allDictAngles[tIdxs[i]].numCodValues;
            if (tScore < minScore && tCNum > 5)
            {
                minScore = tScore;
                rIdx     = i;  
            }
        }
        
    }
    
    void CodClassify::getMatched4dVect(std::vector<std::vector<ID> > tPairs, 
                                       std::vector<ID> aPair, std::vector<int> tIdxs, 
                                       int& rIdx)
    {
        int minScore  = 1000000;
        
        int groupUnit = 100;
        int nonMUnit  = 50;
        int metUnit   = 20;
        
        PeriodicTable tP;
        
        for (int i=0; i < (int)tPairs.size(); i++)
        {
            int tScore = 0;
            tScore += (abs(tP.elements[tPairs[i][0]]["group"]-tP.elements[aPair[0]]["group"])*groupUnit);
            int n1  = StrToInt(tPairs[i][1].substr(3,1));
            int n2  = StrToInt(aPair[1].substr(3,1));
            tScore += (abs(n1-n2)*nonMUnit);
            int m1  = StrToInt(tPairs[i][1].substr(1,1));
            int m2  = StrToInt(aPair[1].substr(1,1));
            tScore += (abs(m1-m2)*metUnit);
            tScore += (abs(tP.elements[tPairs[i][2]]["group"]-tP.elements[aPair[2]]["group"])*groupUnit);
            int n3  = StrToInt(tPairs[i][3].substr(3,1));
            int n4  = StrToInt(aPair[3].substr(3,1));
            tScore += (abs(n3-n4)*nonMUnit);
            int m3  = StrToInt(tPairs[i][3].substr(1,1));
            int m4  = StrToInt(aPair[3].substr(1,1));
            tScore += (abs(m3-m4)*metUnit);
            int tCNum = allDictAngles[tIdxs[i]].numCodValues;
            if (tScore < minScore && tCNum > 5)
            {
                minScore = tScore;
                rIdx     = i;  
            }
        }
        
    }
  
    void CodClassify::getMatched5dVect(std::vector<std::vector<ID> > tPairs, 
                                       std::vector<ID> aPair, std::vector<int> tIdxs, 
                                       int& rIdx)
    {
        int minScore  = 1000000;
        
        int groupUnit = 100;
        int nonMUnit  = 50;
        int metUnit   = 20;
        
        PeriodicTable tP;
        
        for (int i=0; i < (int)tPairs.size(); i++)
        {
            int tScore = 0;
            int n   = StrToInt(tPairs[i][0].substr(3,1));
            int n0  = StrToInt(aPair[0].substr(3,1));
            tScore += (abs(n-n0)*nonMUnit);
            int m   = StrToInt(tPairs[i][0].substr(1,1));
            int m0  = StrToInt(aPair[0].substr(1,1));
            tScore += (abs(m-m0)*metUnit);
            tScore += (abs(tP.elements[tPairs[i][1]]["group"]-tP.elements[aPair[1]]["group"])*groupUnit);
            int n1  = StrToInt(tPairs[i][2].substr(3,1));
            int n2  = StrToInt(aPair[2].substr(3,1));
            tScore += (abs(n1-n2)*nonMUnit);
            int m1  = StrToInt(tPairs[i][2].substr(1,1));
            int m2  = StrToInt(aPair[2].substr(1,1));
            tScore += (abs(m1-m2)*metUnit);
            tScore += (abs(tP.elements[tPairs[i][3]]["group"]-tP.elements[aPair[3]]["group"])*groupUnit);
            int n3  = StrToInt(tPairs[i][4].substr(3,1));
            int n4  = StrToInt(aPair[4].substr(3,1));
            tScore += (abs(n3-n4)*nonMUnit);
            int m3  = StrToInt(tPairs[i][4].substr(1,1));
            int m4  = StrToInt(aPair[4].substr(1,1));
            tScore += (abs(m3-m4)*metUnit);
            int tCNum = allDictAngles[tIdxs[i]].numCodValues;
            if (tScore < minScore && tCNum > 5)
            {
                minScore = tScore;
                rIdx     = i;  
            }
        }
        
    }
    
    void CodClassify::setupTargetAngleUsingdist(std::vector<AngleDict>& tAngles, 
                                                std::vector<AngleDict>::iterator tA, 
                                                int tLev)
    {
        // search for the closest cod-classes of three atoms in an angle
        
        int shortestDist =1000000, iPos=-1;
        
        int dist1=1000000, dist2=1000000;
        
        
        std::vector<std::vector<ID> > targetNBs;
        for (std::vector<int>::iterator iTA=tA->atoms.begin();
                iTA !=tA->atoms.end(); iTA++)
        {
            std::vector<ID> tV;
            StrTokenize(allAtoms[*iTA].codNBSymb, tV, ':');
            targetNBs.push_back(tV);    
        }
        
        // std::cout << (int)targetNBs.size() << std::endl;
        
        // scan the group of angles
        for (int i=0; i <(int)tAngles.size(); i++)
        {
            std::vector<ID> tIdV0, tIdV1, tIdV2;
            
             
            StrTokenize(tAngles[i].atomsNBRep[0], tIdV0, ':');
            StrTokenize(tAngles[i].atomsNBRep[1], tIdV1, ':');
            StrTokenize(tAngles[i].atomsNBRep[2], tIdV2, ':');
            
            
            if ((int)tIdV0.size() == (int)targetNBs[0].size() )
            {
                int dist10 = codAtomsDist(targetNBs[0], tIdV0, tLev);
                
                if ((int)tIdV1.size() == (int)targetNBs[1].size()
                    && (int)tIdV2.size() == (int)targetNBs[2].size())
                {
                    
                    dist1 = dist10 + codAtomsDist(targetNBs[1], tIdV1, tLev) 
                            + codAtomsDist(targetNBs[2], tIdV2, tLev);
                }
                if ((int)tIdV1.size() == (int)targetNBs[2].size()
                    && (int)tIdV2.size() == (int)targetNBs[1].size())
                {
                    dist2 = dist10 + codAtomsDist(targetNBs[1], tIdV2, tLev)
                            + codAtomsDist(targetNBs[2], tIdV1, tLev);
                }
                std::cout << "dist1 " << dist1 << std::endl;  
                std::cout << "dist2 " << dist2 << std::endl; 
            }
            if (dist1 < shortestDist)
            {
                shortestDist = dist1;
                iPos = i;
            }
            if (dist2 < shortestDist)
            {
                shortestDist = dist2;
                iPos = i;
            }
        }
        
        // std::cout << "distance is " << shortestDist << std::endl;
        if (iPos !=-1)
        {
            std::cout << "number of candidate angles " << (int)tAngles.size() << std::endl;
            std::cout << "selected angle position " << iPos << std::endl;
            
            tA->value        = tAngles[iPos].value;
            tA->sigValue     = tAngles[iPos].sigValue;
            if(tA->sigValue > 3.0)
            {
                tA->sigValue = 3.0;
            }
            tA->numCodValues = tAngles[iPos].numCodValues;
            std::cout << "angle value " << tAngles[iPos].value << std::endl;
            std::cout << "Three atom are : " << std::endl;
            for (std::vector<ID>::iterator iA = tAngles[iPos].atomsCodClasses.begin();
                    iA != tAngles[iPos].atomsCodClasses.end(); iA++)
            {
                std::cout << "Atom COD code " 
                          <<  *iA << std::endl;
            }
            
        }
        else
        {
            std::cout << "could not find a similar atom " << std::endl;
        }
                         
       
    }
    
    void CodClassify::setupTargetAngleUsingdist2(std::vector<AngleDict>& tAngles, 
                                                 std::vector<AngleDict>::iterator tA, 
                                                 int as1, int as2, int tLev)
    {
        int a0Idx=0, a1Idx=0, a2Idx=0;
        int iPos = -1, mDiff=1000000;
        std::vector<ID> targetNBI;
        if (tLev==3)
        {
            // the last codNB2Symb
            std::vector<ID> tV2;
            StrTokenize(allAtoms[tA->atoms[as2]].codNB2Symb, tV2, ':'); 
            
            for (int i=0; i < (int)tV2.size(); i++)
            {
                a2Idx+=StrToInt(tV2[i]);
            }
            
            // scan the group of angles
            for (int i=0; i <(int)tAngles.size(); i++)
            {
                std::vector<ID> tAV2;
                int tA2Idx=0; 
                StrTokenize(tAngles[i].atomsNB2Rep[2], tAV2, ':');
                for (int j=0; j < (int)tV2.size(); j++)
                {
                    tA2Idx+=StrToInt(tV2[j]);
                }
                
                int tDiff= abs(tA2Idx-a2Idx);
                
                if (tDiff < mDiff)
                {
                    iPos = i;
                    mDiff = tDiff;
                }
            }            
        }
        else if (tLev==2)
        {
            // Last two NB2 
            // the last codNB2Symb
            std::vector<ID> tV1, tV2;
            StrTokenize(allAtoms[tA->atoms[as1]].codNB2Symb, tV1, ':'); 
            for (int i=0; i < (int)tV1.size(); i++)
            {
                a1Idx+=StrToInt(tV1[i]);
            }
            
            StrTokenize(allAtoms[tA->atoms[as2]].codNB2Symb, tV2, ':'); 
            for (int i=0; i < (int)tV2.size(); i++)
            {
                a2Idx+=StrToInt(tV2[i]);
            }
            
            // scan the group of angles
            for (int i=0; i <(int)tAngles.size(); i++)
            {

                std::vector<ID> tAV1;
                int tA1Idx=0; 
                StrTokenize(tAngles[i].atomsNB2Rep[1], tAV1, ':');
                for (int j=0; j < (int)tV1.size(); j++)
                {
                    tA1Idx+=StrToInt(tV1[j]);
                }
                
                
                std::vector<ID> tAV2;
                int tA2Idx=0; 
                StrTokenize(tAngles[i].atomsNB2Rep[2], tAV2, ':');
                for (int j=0; j < (int)tV2.size(); j++)
                {
                    tA2Idx+=StrToInt(tV2[j]);
                }
                
                int tDiff= abs(tA1Idx-a1Idx) + abs(tA2Idx-a2Idx);
                
                if (tDiff < mDiff)
                {
                    iPos = i;
                    mDiff = tDiff;
                }
            }            
        }
        else if(tLev==1)
        {
            // Last two NB2 
            // the last codNB2Symb
            std::vector<ID> tV0, tV1, tV2;
            StrTokenize(allAtoms[tA->atoms[0]].codNB2Symb, tV0, ':'); 
            for (int j=0; j < (int)tV0.size(); j++)
            {
                a0Idx+=StrToInt(tV0[j]);
            }
            
            StrTokenize(allAtoms[tA->atoms[as1]].codNB2Symb, tV1, ':'); 
            for (int j=0; j < (int)tV1.size(); j++)
            {
                a1Idx+=StrToInt(tV1[j]);
            }
            
            StrTokenize(allAtoms[tA->atoms[as2]].codNB2Symb, tV2, ':'); 
            for (int j=0; j < (int)tV2.size(); j++)
            {
                a2Idx+=StrToInt(tV2[j]);
            }
            
            // scan the group of angles
            for (int i=0; i <(int)tAngles.size(); i++)
            {

                std::vector<ID> tAV0;
                int tA0Idx=0; 
                StrTokenize(tAngles[i].atomsNB2Rep[0], tAV0, ':');
                for (int j=0; j < (int)tV0.size(); j++)
                {
                    tA0Idx+=StrToInt(tV0[j]);
                }
                
                std::vector<ID> tAV1;
                int tA1Idx=0; 
                StrTokenize(tAngles[i].atomsNB2Rep[1], tAV1, ':');
                for (int j=0; j < (int)tV1.size(); j++)
                {
                    tA1Idx+=StrToInt(tV1[j]);
                }
                
                
                std::vector<ID> tAV2;
                int tA2Idx=0; 
                StrTokenize(tAngles[i].atomsNB2Rep[2], tAV2, ':');
                for (int j=0; j < (int)tV2.size(); j++)
                {
                    tA2Idx+=StrToInt(tV2[j]);
                }
                
                int tDiff= abs(tA0Idx-a0Idx)+abs(tA1Idx-a1Idx) + abs(tA2Idx-a2Idx);
                
                if (tDiff < mDiff)
                {
                    iPos = i;
                    mDiff = tDiff;
                }
            }            
        }
           
        
        if (iPos !=-1)
        {
            std::cout << "number of candidate angles " << tAngles.size() << std::endl;
            std::cout << "selected angle position " << iPos << std::endl;
            if (tAngles[iPos].numCodValues >=5)
            {
                tA->value        = tAngles[iPos].value;
                tA->sigValue     = tAngles[iPos].sigValue;
                if(tA->sigValue > 3.0)
                {
                    tA->sigValue = 3.0;
                }
                tA->numCodValues = tAngles[iPos].numCodValues;
                std::cout << "angle value " << tAngles[iPos].value << std::endl;
                std::cout << "Three atom are : " << std::endl;
                for (std::vector<ID>::iterator iA = tAngles[iPos].atomsCodClasses.begin();
                                       iA != tAngles[iPos].atomsCodClasses.end(); iA++)
                {
                    std::cout << "Atom COD code " 
                              <<  *iA << std::endl;
                }
            }
            else
            {
                if (allAtoms[tA->atoms[0]].bondingIdx >0 && 
                    allAtoms[tA->atoms[0]].bondingIdx <4)
                {
                    // std::cout << "Center atom bond index is  " << allAtoms[iAN->atoms[0]].bondingIdx<< std::endl;
                    tA->value = DefaultOrgAngles[allAtoms[tA->atoms[0]].bondingIdx];
                    tA->sigValue = 3.0;
                }
            }
        }
        else
        {
            std::cout << "could not find a similar atom " << std::endl;
        }
        
    }
    
    void CodClassify::setupTargetAngleUsingdistAng(std::vector<AngleDict>& tAngles, 
                                                std::vector<AngleDict>::iterator tA, 
                                                int tLev)
    {
        // search for the closest cod-classes of three atoms in an angle
        
        int shortestDist =1000000, iPos=-1;
        
        int dist1=1000000, dist2=1000000;
        
        
        std::vector<std::vector<ID> > targetNBs;
        for (std::vector<int>::iterator iTA=tA->atoms.begin();
                iTA !=tA->atoms.end(); iTA++)
        {
            std::vector<ID> tV;
            StrTokenize(allAtoms[*iTA].codNBSymb, tV, ':');
            targetNBs.push_back(tV);
        }
        // std::cout << (int)targetNBs.size() << std::endl;
        
        // scan the group of angles
        for (int i=0; i <(int)tAngles.size(); i++)
        {
            std::vector<ID> tIdV0, tIdV1, tIdV2;
           
            
            StrTokenize(tAngles[i].atomsNBRep[0], tIdV0, ':');
            StrTokenize(tAngles[i].atomsNBRep[1], tIdV1, ':');
            StrTokenize(tAngles[i].atomsNBRep[2], tIdV2, ':');
            
            
            if ((int)tIdV0.size() == (int)targetNBs[0].size() )
            {
                int dist10 = codAtomsDist(targetNBs[0], tIdV0, tLev);
                
                if ((int)tIdV1.size() == (int)targetNBs[1].size()
                    && (int)tIdV2.size() == (int)targetNBs[2].size())
                {
                    
                    dist1 = dist10 + codAtomsDist(targetNBs[1], tIdV1, tLev) 
                            + codAtomsDist(targetNBs[2], tIdV2, tLev);
                }
                if ((int)tIdV1.size() == (int)targetNBs[2].size()
                    && (int)tIdV2.size() == (int)targetNBs[1].size())
                {
                    dist2 = dist10 + codAtomsDist(targetNBs[1], tIdV2, tLev)
                            + codAtomsDist(targetNBs[2], tIdV1, tLev);
                }
                // std::cout << "dist1 " << dist1 << std::endl;  
                // std::cout << "dist2 " << dist2 << std::endl; 
            }
            if (dist1 < shortestDist)
            {
                shortestDist = dist1;
                iPos = i;
            }
            if (dist2 < shortestDist)
            {
                shortestDist = dist2;
                iPos = i;
            }
        }
        
        // std::cout << "distance is " << shortestDist << std::endl;
        if (iPos !=-1)
        {
            tA->value        = tAngles[iPos].value;
            tA->sigValue     = tAngles[iPos].sigValue;
            if(tA->sigValue > 3.0)
            {
                tA->sigValue = 3.0;
            }
            
            tA->numCodValues = tAngles[iPos].numCodValues;
            // std::cout << "angle value " << tAngles[iPos].value << std::endl;
            // std::cout << "Three atom are : " << std::endl;
            //for (std::vector<int>::iterator iA = tAngles[iPos].atoms.begin();
            //        iA != tAngles[iPos].atoms.end(); iA++)
            // {
            //    std::cout << "Atom " << allAtoms[*iA].id << " COD code " 
            //              <<  allAtoms[*iA].codClass << std::endl;
            // }
            
        }
        else
        {
            std::cout << "could not find a similar atom " << std::endl;
        }                    
       
    }
    
    void CodClassify::setupTargetAngleUsingMean(std::vector<AngleDict>& tAngles, 
                                                std::vector<AngleDict>::iterator tAN)
    {
        REAL tValueSum=0.0, tSqValueSum=0.0,  tMeanSqValue=0.0;
        REAL tMeanValue=0.0, tSqMeanValue=0.0, tStVar=0.0;
        int  nTot=0;
        
        for (std::vector<AngleDict>::iterator iAN=tAngles.begin();
                iAN !=tAngles.end(); iAN++)
        {
            nTot +=iAN->numCodValues;
            tValueSum   += (iAN->value*(iAN->numCodValues));
            // std::cout << "a angle value " << iAN->value << std::endl;
            tSqValueSum += (iAN->value*(iAN->value)*(iAN->numCodValues));
        }
        
        
        tMeanValue = tValueSum/nTot;
        tSqMeanValue = (tMeanValue*tMeanValue);
        // std::cout << "Mean " <<tMeanValue <<std::endl;
        // std::cout << "Number of values for average " << nTot << std::endl;
        tMeanSqValue  = (tSqValueSum/nTot);
        
        tStVar  = std::fabs(tMeanSqValue-tSqMeanValue);
        tStVar  = std::sqrt(tStVar);
        
        tAN->value    = tMeanValue;
        tAN->sigValue = tStVar;
        if (tAN->sigValue > 3.0)
        {
            tAN->sigValue = 3.0;
        }
        tAN->numCodValues = nTot;

    }
    
    void CodClassify::setupTargetAngleUsingMean2(std::vector<AngleDict>& tAngles, 
                                                std::vector<AngleDict>::iterator tAN,
                                                int tHa1, int tHa2, int tHa3,
                                                ID tA1NB2, ID tA2NB2, ID tA3NB2)
    {
        REAL tValueSum=0.0, tSqValueSum=0.0,  tMeanSqValue=0.0;
        REAL tMeanValue=0.0, tSqMeanValue=0.0, tStVar=0.0;
        int  nTot=0;
        
        for (std::vector<AngleDict>::iterator iAN=tAngles.begin();
                iAN !=tAngles.end(); iAN++)
        {
            nTot +=iAN->numCodValues;
            tValueSum   += (iAN->value*(iAN->numCodValues));
            // std::cout << "a angle value " << iAN->value << std::endl;
            tSqValueSum += (iAN->value*(iAN->value)*(iAN->numCodValues));
        }
        
        if (nTot >=10)
        {
            tMeanValue = tValueSum/nTot;
            tSqMeanValue = (tMeanValue*tMeanValue);
            // std::cout << "Mean " <<tMeanValue <<std::endl;
            // std::cout << "Number of values for average " << nTot << std::endl;
            tMeanSqValue  = (tSqValueSum/nTot);
        
            tStVar  = std::fabs(tMeanSqValue-tSqMeanValue);
            tStVar  = std::sqrt(tStVar);
        
            tAN->value    = tMeanValue;
            tAN->sigValue = tStVar;
            if (tAN->sigValue > 3.0)
            {
                tAN->sigValue = 3.0;
            }
            tAN->numCodValues = nTot;
        }
        else
        {
            tAN->value        = allDictAnglesIdx2[tHa1][tHa2][tHa3][tA1NB2][tA2NB2][tA3NB2][0].value;
            tAN->sigValue     = allDictAnglesIdx2[tHa1][tHa2][tHa3][tA1NB2][tA2NB2][tA3NB2][0].sigValue;
            tAN->numCodValues = allDictAnglesIdx2[tHa1][tHa2][tHa3][tA1NB2][tA2NB2][tA3NB2][0].numCodValues;
        }
    }
   
    void CodClassify::setupTargetAngleUsingMean3(std::vector<REAL>  & tAngValues,
                                        std::vector<int>            & tAngNums,
                                        std::vector<AngleDict>::iterator tAN)
                                                 
    {
        REAL tValueSum=0.0, tSqValueSum=0.0,  tMeanSqValue=0.0;
        REAL tMeanValue=0.0, tSqMeanValue=0.0, tStVar=0.0;
        int  nTot=0;
        
        for (unsigned i=0; i < tAngValues.size(); i++)
        {
            nTot+=tAngNums[i];
            REAL tV = tAngValues[i]*tAngNums[i];
            tValueSum   += tV;
            tSqValueSum += (tAngValues[i]*tV);
        }
        
        if (nTot >=1)
        {
            tMeanValue = tValueSum/nTot;
            tSqMeanValue = (tMeanValue*tMeanValue);
            // std::cout << "Mean " <<tMeanValue <<std::endl;
            // std::cout << "Number of values for average " << nTot << std::endl;
            tMeanSqValue  = (tSqValueSum/nTot);
        
            tStVar  = std::fabs(tMeanSqValue-tSqMeanValue);
            tStVar  = std::sqrt(tStVar);
        
            tAN->value    = tMeanValue;
            tAN->sigValue = tStVar;
            if (tAN->sigValue > 3.0)
            {
                tAN->sigValue = 3.0;
            }
            tAN->numCodValues = nTot;
        }
        else
        {
            std::cout << "Number of angle values involved in mean calculation is zero "
                      << std::endl;
            exit(1);
        }
    }
    
    void CodClassify::setupTargetAngles()
    {
        //initTargetAngles();
        //std::cout << "Finished initial target angles " << std::endl;  
        setDefaultOrgAngle();
        setDefaultCoordGeos();
        std::cout << "Finish setDefaultCoordGeos() " << std::endl;
        groupCodOrgAngles2();
        std::cout << "Finish groupCodOrgAngles() " << std::endl;
        groupCodMetAngles();
        std::cout << "Finish groupCodMetAngles() " << std::endl;
        
        groupCodAnglesWithNonCenteredMetal();
        std::cout << "Finish groupCodAnglesWithNonCenteredMetal() " << std::endl;
        searchCodAngles();
        
        std::cout << "Finish searching angles " << std::endl;
        /*
        for (std::vector<AngleDict>::iterator iA=allAngles.begin();
                iA!=allAngles.end(); iA++)
        {
            std::cout << "For angle formed by atom " << allAtoms[iA->atoms[0]].id << "(Center)"
                      << " atom " << allAtoms[iA->atoms[1]].id 
                      << " atom " << allAtoms[iA->atoms[2]].id << " : " << std::endl
                      << "Value " << iA->value << std::endl
                      << "sigValue " << iA->sigValue << std::endl;  
        }
         */
        
    }
    
    void CodClassify::setupTargetAngles2()
    {
        //initTargetAngles();
        //std::cout << "Finished initial target angles " << std::endl;  
        setDefaultOrgAngle();
        setDefaultCoordGeos();
        std::cout << "Finish setDefaultCoordGeos() " << std::endl;
        groupCodOrgAngles22();
        std::cout << "Finish groupCodOrgAngles() " << std::endl;
        groupCodMetAngles();
        std::cout << "Finish groupCodMetAngles() " << std::endl;
        
        groupCodAnglesWithNonCenteredMetal();
        std::cout << "Finish groupCodAnglesWithNonCenteredMetal() " << std::endl;
        searchCodAngles2();
        
        std::cout << "Finish searching angles " << std::endl;
        /*
        for (std::vector<AngleDict>::iterator iA=allAngles.begin();
                iA!=allAngles.end(); iA++)
        {
            std::cout << "For angle formed by atom " << allAtoms[iA->atoms[0]].id << "(Center)"
                      << " atom " << allAtoms[iA->atoms[1]].id 
                      << " atom " << allAtoms[iA->atoms[2]].id << " : " << std::endl
                      << "Value " << iA->value << std::endl
                      << "sigValue " << iA->sigValue << std::endl;  
        }
         */
        
    }
    /*
    void CodClassify::setupTargetAnglesUsingSqlite()
    {
                //initTargetAngles();
        //std::cout << "Finished initial target angles " << std::endl;  
        setDefaultOrgAngle();
        setDefaultCoordGeos();
        std::cout << "Finish setDefaultCoordGeos() " << std::endl;
        // groupCodOrgAngles();
        // std::cout << "Finish groupCodOrgAngles() " << std::endl;
        groupCodMetAngles();
        std::cout << "Finish groupCodMetAngles() " << std::endl;
        groupCodAnglesWithNonCenteredMetal();
        std::cout << "Finish groupCodAnglesWithNonCenteredMetal() " << std::endl;
        // searchCodAngles();
        searchCodAnglesUsingSqlite();
        std::cout << "Finish searching angles " << std::endl;
        
    }
   */
    
    // Torsion angles related 
    
    bool CodClassify::checkATorsAtomsInPla(std::vector<int> & tAtms)
    {
        
        for (std::vector<PlaneDict>::iterator iP=allPlanes.begin();
                        iP !=allPlanes.end(); iP++)
        {
            std::vector<ID> inPlAtoms;
            for (std::vector<int>::iterator iA=tAtms.begin();
                    iA !=tAtms.end(); iA++)
            {
                std::map<ID, int>::iterator iFind =iP->atoms.find(allAtoms[*iA].id);
                if (iFind !=iP->atoms.end())
                {
                    inPlAtoms.push_back(allAtoms[*iA].id);
                }
            }
            
            if ((int)inPlAtoms.size() == (int)tAtms.size())
            {
                return true;
            }
        }
        
        
        return false;
    }
    
    void CodClassify::fixTorIDs()
    {
        int idxTors  = 1, idxPTors = 1, idxSp3Sp3=1, idxSp2Sp3=1, idxSp2Sp2=1;
        //std::cout << "There are " << (int)allTorsions.size() 
        //          << "Torsions" << std::endl;
       
        
        for (std::vector<TorsionDict>::iterator iT=allTorsions.begin();
                        iT !=allTorsions.end(); iT++)
        {
            // std::cout << "look at torsion " << iT->seriNum << std::endl;
            if(checkATorsAtomsInPla(iT->atoms))
            {
                iT->id = "const_sp2_sp2_" + IntToStr(idxPTors);
                if (iT->id.size() >=16 )
                {
                    iT->id = "const_" + IntToStr(idxPTors);
                }
                //iT->id = "P_sp2_sp2_" + IntToStr(idxPTors);
                idxPTors +=1;
            }
            else
            {
                if (iT->period ==3 )
                {
                    iT->id = "sp3_sp3_" + IntToStr(idxSp3Sp3);
                    idxSp3Sp3+=1;
                }
                else if (iT->period == 6)
                {
                    iT->id= "sp2_sp3_"+IntToStr(idxSp2Sp3);
                    idxSp2Sp3+=1;
                }
                else if (iT->period == 2)
                {
                    iT->id = "sp2_sp2_"+IntToStr(idxSp2Sp2);
                    idxSp2Sp2+=1;
                    
                }
                else
                {
                    iT->id = "other_tor_" + IntToStr(idxTors);
                    idxTors++;
                }
            }
            //std::cout << "its ID now is " << iT->id << std::endl; 
        }
        
        
    }
    
    void CodClassify::checkAngConstraints()
    {
        std::map<int, std::vector<int> >  tSP2Angs, tSP3Angs;
        
        
        std::cout << "There are " << allRings.size() << " rings " << std::endl
                  << "These are: " << std::endl;
        for (std::map<ID, std::vector<RingDict> >::iterator iR=allRings.begin();
                iR !=allRings.end(); iR++)
        {
            std::cout << "Ring " << iR->first << "." << std::endl;
            for (std::vector<RingDict>::iterator iRV=iR->second.begin();
                    iRV !=iR->second.end(); iRV++)
            {
                std::cout << "Ring size:  " << iRV->atoms.size() << std::endl;;
                for (std::vector<AtomDict>::iterator iA=iRV->atoms.begin(); 
                        iA != iRV->atoms.end(); iA++)
                {
                    std::cout<< "Atom " << iA->id << std::endl;
                }
                checkRingAngleConstraints(iRV);
            }
        }
        
        
        // Check all sp2 and sp3 angles centered at different atoms
        
        for (unsigned i=0;  i < allAngles.size(); i++)
        {
            int iCenAtom= allAngles[i].atoms[0];
            if (allAtoms[iCenAtom].bondingIdx==2)
            {
                tSP2Angs[iCenAtom].push_back(i);
            }
            else if (allAtoms[iCenAtom].bondingIdx==3)
            {
                tSP3Angs[iCenAtom].push_back(i);
            }
        }
        
        // check all sp2 atoms
        for (std::map<int, std::vector<int> >::iterator iSetAngs=tSP2Angs.begin();
                iSetAngs !=tSP2Angs.end(); iSetAngs++)
        {
            checkSP2Constraints(iSetAngs->second);
        }
        
        // check all sp3 atoms
        for (std::map<int, std::vector<int> >::iterator iSetAngs=tSP3Angs.begin();
                iSetAngs !=tSP3Angs.end(); iSetAngs++)
        {
            checkSP3Constraints(iSetAngs->second);
        }
        
    }
    
    void CodClassify::checkRingAngleConstraints(std::vector<RingDict>::iterator tRv)
    {
        std::vector<ID> atomIDs;
        std::vector<int> angIdxs;
        std::map<int, REAL> moF;
        moF[0] = 0.0;
        moF[1] = 0.25;
        moF[2] = 0.50;
        moF[3] = 0.75;
        for (int i=4; i <=10; i++)
        {
            moF[i] = 1.0;
        }
        
        if (tRv->atoms.size()==0)
        {
            std::cout << "Ring " << tRv->rep << " contains no atoms , Bug?"
                    << std::endl;
            return;
        }
        
        
        for (std::vector<AtomDict>::iterator iA=tRv->atoms.begin();
                iA!=tRv->atoms.end(); iA++)
        {
            if (iA->bondingIdx == 3)
            {
                // Non-plane ring, 
                return;
            }
            
            atomIDs.push_back(iA->id);
        }
        
        for (unsigned i=0; i < allAngles.size(); i++)
        {
            if (std::find(atomIDs.begin(), atomIDs.end(), allAtoms[allAngles[i].atoms[0]].id)
                    !=atomIDs.end())
            {
                if(std::find(atomIDs.begin(), atomIDs.end(), allAtoms[allAngles[i].atoms[1]].id)
                    !=atomIDs.end() && 
                    std::find(atomIDs.begin(), atomIDs.end(), allAtoms[allAngles[i].atoms[2]].id)
                    !=atomIDs.end())
                {
                    angIdxs.push_back(i);
                    /*
                    std::cout << "Initial angle: center atom " 
                              <<  allAtoms[allAngles[i].atoms[0]].id
                              << " atom1 " << allAtoms[allAngles[i].atoms[1]].id
                              << " atom2 " << allAtoms[allAngles[i].atoms[2]].id
                              << " value " << allAngles[i].value 
                              << " sigma " << allAngles[i].sigValue << std::endl;
                     */
                }
                
            }
        }
        
        // check 
        if (angIdxs.size() !=atomIDs.size())
        {
            std::cout << "Some of ring angles are missing, Bug? " << std::endl;
            exit(1);
        }
        else
        {
            int rSize = (int)atomIDs.size();
            REAL aSum=0.0, aM=(rSize-2)*180.0/rSize, aMDev=0.0;
            int aSize =0;
            REAL curDev = 0.0;
            // First round 
            for (std::vector<int>::iterator iAd=angIdxs.begin();
                    iAd !=angIdxs.end(); iAd++)
            {
                
                if (allAngles[*iAd].levelCodValue==0)
                {
                    allAngles[*iAd].isFixed = true;
                }
                else
                {
                    aSum+=allAngles[*iAd].value;
                    aSize++;
                }
            }
            
            if (aSize !=0)
            {
                aMDev = aM-aSum/aSize;
            
                for (std::vector<int>::iterator iIdx=angIdxs.begin();
                        iIdx !=angIdxs.end(); iIdx++)
                {
                    if(!allAngles[*iIdx].isFixed)
                    {
                        curDev=allAngles[*iIdx].value-aM;
                        allAngles[*iIdx].value-=(moF[allAngles[*iIdx].levelCodValue]*curDev);
                            
                       
                    }
                    
                }
            }
            // Second round 
            
            
            aSum=0.0;
            for (std::vector<int>::iterator iAd=angIdxs.begin();
                    iAd !=angIdxs.end(); iAd++)
            {
                aSum+=allAngles[*iAd].value;
            }
            
            aMDev = aM -aSum/rSize;
            
            for (std::vector<int>::iterator iAd=angIdxs.begin();
                    iAd !=angIdxs.end(); iAd++)
            {
                allAngles[*iAd].value +=(aMDev);
                allAngles[*iAd].isFixed = true;
            }
        }
        
        // output the angle values
       
        // std::cout << "Internal angles in the ring " << tRv->rep << std::endl;
        
        REAL aSum=0.0;
        for (std::vector<int>::iterator iAd=angIdxs.begin();
               iAd !=angIdxs.end(); iAd++)
        {
            aSum+=allAngles[*iAd].value;
            std::cout << "Angle between " << allAtoms[allAngles[*iAd].atoms[0]].id 
                      << "(center) and " << allAtoms[allAngles[*iAd].atoms[1]].id
                      << " and " << allAtoms[allAngles[*iAd].atoms[2]].id
                      << " value: " << allAngles[*iAd].value << std::endl;
                    
        }
        std::cout << "Sum for the values of all internal angles " 
                  << aSum << std::endl;
       
    }
    
    void CodClassify::checkSP2Constraints(std::vector<int> tAngIdxs)
    {
        
        if (tAngIdxs.size()==3)
        {
            REAL sp2Sum;
            REAL sp2Diff =0.0;
            REAL tAngSum = 0.0;
            int  nNonF=0;
            REAL vF   =0.0;
            for (std::vector<int>::iterator iA=tAngIdxs.begin();
                    iA !=tAngIdxs.end(); iA++)
            {
                if (!allAngles[*iA].isFixed)
                {
                    tAngSum+=(allAngles[*iA].value);
                    nNonF+=1;
                }
                else
                {
                    vF+=(allAngles[*iA].value);
                }
            }
            
            if (nNonF >0)
            {
                sp2Sum  = (360.0-vF);
                sp2Diff = (sp2Sum - tAngSum)/nNonF;
            
                if (fabs(sp2Diff) >0.01)
                {
                    // add the diff to individual angles
                    tAngSum=0.0;
                    for (std::vector<int>::iterator iA=tAngIdxs.begin();
                            iA !=tAngIdxs.end(); iA++)
                    {
                        if (!allAngles[*iA].isFixed)
                        {
                            allAngles[*iA].value +=sp2Diff;
                            tAngSum+=(allAngles[*iA].value);
                        }
                    }
                    // Further check.  
                    sp2Diff = sp2Sum - tAngSum;
                    //  Now sp2Diff should be negligible small, add it to an angle
                    allAngles[tAngIdxs[0]].value += sp2Diff;
                }
                else 
                {
                    // very small, just add it to an angle 
                    allAngles[tAngIdxs[0]].value += sp2Diff;
                }
            }
            
            //Debug, output  the modified sum, should be 360.00
            /*
            std::cout << "values for angles centered at atom "
                      << allAtoms[allAngles[tAngIdxs[0]].atoms[0]].id 
                      << std::endl;
            tAngSum=0.0;
            for (std::vector<int>::iterator iA=tAngIdxs.begin();
                        iA !=tAngIdxs.end(); iA++)
            {
                std::cout << "Angle between atoms " 
                          << allAtoms[allAngles[*iA].atoms[0]].id
                          << "(center), " << allAtoms[allAngles[*iA].atoms[1]].id
                          << " and " << allAtoms[allAngles[*iA].atoms[2]].id
                          << ", value: "<< allAngles[*iA].value << std::endl;
                tAngSum+=(allAngles[*iA].value);
            }
            
            std::cout << "Sum of values for angles centered at atom " 
                      << allAtoms[allAngles[tAngIdxs[0]].atoms[0]].id 
                      << " is " << tAngSum << std::endl;
            */
        }
    }
    
    void CodClassify::checkSP3Constraints(std::vector<int> tAngIdxs)
    {
    }
    
    bool CodClassify::checkSpeAng(std::vector<AngleDict>::iterator tA)
    {
        // at the moment, check if the center
        // atom of the angle is B with 5 connections
        // or C with 6 connections 
        if ((allAtoms[tA->atoms[0]].chemType == "C" && allAtoms[tA->atoms[0]].connAtoms.size()==6)
             || (allAtoms[tA->atoms[0]].chemType == "B" && allAtoms[tA->atoms[0]].connAtoms.size()==5))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    
    void CodClassify::setSpecialAngles(std::map<int, std::vector<AngleDict> > & tAngs)
    {
        // If atom is B with 5 connections or C with 6 connections,
        // we do not search the database but set the values of those angles
        // via special configuration the atoms are in
        
        for (std::map<int, std::vector<AngleDict> >::iterator iAs=tAngs.begin();
                iAs !=tAngs.end(); iAs++)
        {
            if (allAtoms[iAs->first].chemType=="B" 
                && allAtoms[iAs->first].connAtoms.size()==5)
            {
                setOneSetBoronAngles(iAs);
            }
            else if (allAtoms[iAs->first].chemType=="C" 
                     && allAtoms[iAs->first].connAtoms.size()==6)
            {
                setOneSetCarbonAngles(iAs);
            }
            else
            {
                std::cout << "Bug:: This set of angles should not appear in function setSpecialAngles() "
                          << std::endl << "This set of angles centered about atom "
                          << allAtoms[iAs->first].id << ". It bonds to " 
                          << allAtoms[iAs->first].connAtoms.size() << " atoms "
                          << std::endl;
                exit(1);
            }
        }
    }
    
    void CodClassify::setOneSetBoronAngles(std::map<int,std::vector<AngleDict> >::iterator tAs)
    {
        // Deal with Boron with 5 bonds
        for (std::vector<AngleDict>::iterator iA=tAs->second.begin();
                iA !=tAs->second.end(); iA++)
        {
            if (std::find(allAtoms[iA->atoms[1]].connAtoms.begin(), 
                          allAtoms[iA->atoms[1]].connAtoms.end(), iA->atoms[2]) 
                    != allAtoms[iA->atoms[1]].connAtoms.end())
            {
                // two atoms connected within a triangle
                iA->value = 60.0;
                iA->valueST = iA->value;
                iA->sigValue = 3.0;
            }
            else
            {
                // Two are not connected. There is only option
                iA->value = 108.0;
                iA->valueST = iA->value;
                iA->sigValue = 3.0;
            }
        }
    }
    
    void CodClassify::setOneSetCarbonAngles(std::map<int,std::vector<AngleDict> >::iterator tAs)
    {
        // Deal with Carbon atoms with 6 bonds
        // 1. In all atoms the center atom connected, there should be one atom
        //    that does not any other atoms (except the center atom). Find that one
        int aUp=-1;
        for (unsigned i=0; i < allAtoms[tAs->first].connAtoms.size(); i++)
        {
            aUp =allAtoms[tAs->first].connAtoms[i];
            for (unsigned j=i+1; j < allAtoms[tAs->first].connAtoms.size(); j++)
            {
                if (std::find(allAtoms[allAtoms[tAs->first].connAtoms[j]].connAtoms.begin(), 
                              allAtoms[allAtoms[tAs->first].connAtoms[j]].connAtoms.end(),
                              aUp) !=allAtoms[allAtoms[tAs->first].connAtoms[j]].connAtoms.end())
                {
                    // aUp is connected to allAtoms[tAs->first].connAtoms[j]
                    aUp = -1;
                    break;
                }
            }
            if (aUp !=-1)
            {
                break;
            }
        }
        
        if (aUp !=-1)
        {
            // These values are for top angles
            REAL A1 = 54.0 * PI180;
            REAL C1 = 1.0/(2.0*cos(A1));
            REAL T1 = 90.0 + acos(C1)*PID180;
            REAL dT = 3.0;
            for (std::vector<AngleDict>::iterator iA=tAs->second.begin();
                iA !=tAs->second.end(); iA++)
            {
                if (std::find(iA->atoms.begin(), iA->atoms.end(), aUp)
                    ==iA->atoms.end())
                {
                    // The top atom is not at this angle 
                    if (std::find(allAtoms[iA->atoms[1]].connAtoms.begin(), 
                        allAtoms[iA->atoms[1]].connAtoms.end(), iA->atoms[2]) 
                        != allAtoms[iA->atoms[1]].connAtoms.end())
                    {
                        // Two atoms connected within a triangle
                        iA->value    = 60.0;
                        iA->valueST  = iA->value;
                        iA->sigValue = dT; 
                    }
                    else
                    {
                        // Two are not connected. There is only option
                        iA->value    = 108.0;
                        iA->valueST  = iA->value;
                        iA->sigValue = dT;
                    }
                }
                else
                {
                    //This angle contain the top atom
                    iA->value    = T1;
                    iA->valueST  = iA->value;
                    iA->sigValue = dT;
                }
            }
        }
        
    }
    void CodClassify::setupTargetTorsions()
    {
        // search COD for torsion angle values, no needed at the moments.
    }
    
    // All target values related  
    
    void CodClassify::setupAllTargetValues()
    {
        
        
        // std::cout << "libmol table should be " << libmolTabDir << std::endl;
       
        setupTargetBonds();
        // setupTargetBondsUsingSqlite();
        
        setupTargetAngles();
        // setupTargetAnglesUsingSqlite();
        // Torsion angle values have been setup when atoms were read in
        // from an input cif file. The following is in case we need the 
        // values from COD
        // setupTargetTorsions();
       
        
        fixTorIDs();
        
        // Chiral centers have been setup when atoms are read in from an
        // a input cif file
        
        // Initiate a set of rough coordinates of atoms
        // initRoughCoords();
        setupAllStdValues();
        
    }
    
    void CodClassify::setupAllTargetValues2()
    {
        
        
        // std::cout << "libmol table should be " << libmolTabDir << std::endl;
        
        setupTargetBonds2();
        // setupTargetBondsUsingSqlite();
        
        setupTargetAngles2();
        // setupTargetAnglesUsingSqlite();
        // Torsion angle values have been setup when atoms were read in
        // from an input cif file. The following is in case we need the 
        // values from COD
        // setupTargetTorsions();
        
        fixTorIDs();
        
        // Chiral centers have been setup when atoms are read in from an
        // a input cif file
        
        // Initiate a set of rough coordinates of atoms
        // initRoughCoords();
        setupAllStdValues();
        
    }
    
    
    
    
    void CodClassify::setupAllStdValues()
    {
        //std::cout << "Set all bond std values " << std::endl; 
        for (std::vector<BondDict>::iterator iB=allBonds.begin();
                iB != allBonds.end(); iB++)
        {
            iB->valueST = iB->value;
            if (iB->sigValue < 0.01 || iB->sigValue > 0.02)
            {
                iB->sigValue= 0.02;
            }
            iB->sigValueST = iB->sigValue;
            
        }
        //std::cout << "Set all angle std values " << std::endl; 
        for (std::vector<AngleDict>::iterator iAn = allAngles.begin();
                iAn !=allAngles.end(); iAn++)
        {
            iAn->valueST    = iAn->value;
            if (iAn->sigValue < 0.1 || iAn->sigValue > 3.0)
            {
                iAn->sigValue = 3.0;
            }
            iAn->sigValueST = iAn->sigValue;
        }
        
        //std::cout << "Set all torsion std values " << std::endl; 
        for (std::vector<TorsionDict>::iterator iTo = allTorsions.begin();
                iTo != allTorsions.end(); iTo++)
        {
            iTo->valueST    = iTo->value;
            iTo->sigValueST = iTo->sigValue;
        }
        
        //std::cout << "Set all chiral std values " << std::endl; 
        for (std::vector<ChiralDict>::iterator iCh=allChirals.begin();
                iCh !=allChirals.end(); iCh++)
        {
            iCh->valueST = iCh->value;
            iCh->signST  = iCh->sign;
        }
        
        
        
    }
    
    void CodClassify::initRoughCoords()
    {
        
        std::cout << "Initiate atomic coordinates " << std::endl;
        initRoughCoordsTor();
        
        int dim = (int)allAtoms.size();
       
        REAL tDistMat[dim][dim];
        
        for (int i=0; i < dim; i++)
        {
            for (int j=0; j < dim; j++)
            {
                tDistMat[i][j] =0.0;
            }
        }
        
        // initiate the elements of the distance matrix by bonds
        for (std::vector<BondDict>::iterator iB=allBonds.begin();
                iB !=allBonds.end(); iB++)
        {
            std::vector<int> tPos;
            for (std::map<ID, int>::iterator iM=iB->fullAtoms.begin();
                    iM != iB->fullAtoms.end(); iM++)
            {
                tPos.push_back(iM->second);
            }
            
            tDistMat[tPos[0]][tPos[1]] = iB->value;
            tDistMat[tPos[1]][tPos[0]] = iB->value;
        }
        
        std::map<ID, std::vector<AngleDict> > metAngles;
        // initiate the elements of the distance matrix contributed by angles 
        
        // Dealt with organic element centered angles and store 
        // metal centered angles in a vector
        for (std::vector<AngleDict>::iterator iA=allAngles.begin();
                iA != allAngles.end(); iA++)
        {
            int i1, i2, i3;
            
            std::sort(iA->atoms.begin()+1, iA->atoms.end());
            i1 = iA->atoms[0]; // the center atom of the angle
            i2 = iA->atoms[1];
            i3 = iA->atoms[2];
            
            if (allAtoms[i1].isMetal)
            {    
                ID tID = IntToStr(i2) + "_" + IntToStr(i1) + "_" + IntToStr(i3);
                metAngles[tID].push_back(*iA);
            }
            else
            {
                REAL dist23 = sqrt(tDistMat[i2][i1]*tDistMat[i2][i1] 
                              +tDistMat[i1][i3]*tDistMat[i1][i3]
                              -2.0*tDistMat[i2][i1]*tDistMat[i1][i3]*cos(iA->value*PI180));
            
               if (fabs(tDistMat[i2][i3]) <=0.001) // prevent the angle values mixed 
               {
                   tDistMat[i2][i3]= dist23;
                   tDistMat[i3][i2]= dist23;
                
                   /*
                   std::cout << "Angle between " << i2 << "-" << i1 << "-" << i3 
                             << " is " << iA->value << " in degrees and "
                             << iA->value*PI180 << " in rads " << std::endl;
                   std::cout << "dist between " << allAtoms[i2].id  << " and "
                             << allAtoms[i3].id  << " is " << tDistMat[i2][i3] << std::endl;
                   */
               }
            }
        }
        
        // Now, dealt with multiple values of each metal centered angles
        for (std::map<ID, std::vector<AngleDict> >::iterator iMA=metAngles.begin();
                iMA !=metAngles.end(); iMA++)
        {
            int i1 = iMA->second[0].atoms[0]; // the center atom of the angle
            int i2 = iMA->second[0].atoms[1];
            int i3 = iMA->second[0].atoms[2];
            
            if ((int)iMA->second.size() >1)
            {
                // multiple values of the same angle
            
                REAL dist23_1 = sqrt(tDistMat[i2][i1]*tDistMat[i2][i1] 
                                +tDistMat[i1][i3]*tDistMat[i1][i3]
                                -2.0*tDistMat[i2][i1]*tDistMat[i1][i3]*cos(iMA->second[0].value*PI180));
                REAL dist23_2 = sqrt(tDistMat[i2][i1]*tDistMat[i2][i1] 
                                +tDistMat[i1][i3]*tDistMat[i1][i3]
                                -2.0*tDistMat[i2][i1]*tDistMat[i1][i3]*cos(iMA->second[1].value*PI180));
               if (dist23_1 >=dist23_2)
               {
                   tDistMat[i2][i3]= dist23_1;
                   tDistMat[i3][i2]= dist23_2;
               }
               else
               {
                   tDistMat[i2][i3]= dist23_2;
                   tDistMat[i3][i2]= dist23_1;
               }
            }
            else
            {
                REAL dist23 = sqrt(tDistMat[i2][i1]*tDistMat[i2][i1] 
                               +tDistMat[i1][i3]*tDistMat[i1][i3]
                               -2.0*tDistMat[i2][i1]*tDistMat[i1][i3]*cos(iMA->second[0].value*PI180));
                if (fabs(tDistMat[i2][i3]) <=0.001) // prevent the angle values mixed 
                {
                   tDistMat[i2][i3]= dist23;
                   tDistMat[i3][i2]= dist23;
                }
            }
        }
        
        // initiate the elements of distance matrix contributed by torsion 
        
        REAL cosmin = 1.0;
        REAL cosmax = -1.0;
        for (std::vector<TorsionDict>::iterator iTo=allTorsions.begin();
                iTo !=allTorsions.end(); iTo++)
        {
            int i1 = iTo->atoms[0]; // the center atom of the angle
            int i2 = iTo->atoms[1];
            int i3 = iTo->atoms[2];
            int i4 = iTo->atoms[3];
            
            if (tDistMat[i1][i4] < 0.0001)
            {
                REAL rab,rbc,rac,rcd,rbd;
                REAL cosabc, sinabc, cosbcd, sinbcd, bndmin, bndmax;             
                rab = tDistMat[i2][i1];
                rbc = tDistMat[i3][i2];
                rcd = tDistMat[i4][i3];
            
                rac = tDistMat[i3][i1];
                rbd = tDistMat[i4][i2];
                if (rab > 0.0001 && rbc > 0.0001 && rcd > 0.0001)
                {
                    cosabc = (rab*rab + rbc*rbc - rac*rac)/(2.0*rab*rbc);
                    sinabc = sqrt(std::max(0.0,(1.0 - cosabc*cosabc)));
                    cosbcd = (rbc*rbc + rcd*rcd - rbd*rbd)/(2.0*rbc*rcd);
                    sinbcd = sqrt(std::max(0.0,(1.0 - cosbcd*cosbcd)));
                    bndmin = rab*rab + rbc*rbc + rcd*rcd + 2.0*rab*rcd*cosabc*cosbcd
                             - 2.0*rab*rcd*sinabc*sinbcd*cosmin
                             - 2.0*rbc*(rab*cosabc+rcd*cosbcd);
                    bndmax = rab*rab + rbc*rbc + rcd*rcd + 2.0*rab*rcd*cosabc*cosbcd
                             - 2.0*rab*rcd*sinabc*sinbcd*cosmax
                             - 2.0*rbc*(rab*cosabc+rcd*cosbcd);
                    bndmin = sqrt(std::max(0.0,bndmin));
                    bndmax = sqrt(std::max(0.0,bndmax));
                    
                    if (i1 > i4)
                    {
                        tDistMat[i1][i4] = bndmin;
                        tDistMat[i4][i1] = bndmax;
                    }
                    else
                    {
                        tDistMat[i1][i4] = bndmax;
                        tDistMat[i4][i1] = bndmin;
                    }
                }          
            } 
        }
        /*
        // Set the rest of matrix to a long distance
        for (int i=0; i < dim; i++)
        {
            for (int j=i+1; j < dim; j++)
            {
                if (fabs(tDistMat[i][j]) < 0.0001)
                {
                    tDistMat[i][j] =(int)allAtoms.size()*1.5;
                    tDistMat[j][i] =2.5;
                }
            }
        }
        */
        // initiate random coordinates 
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA !=allAtoms.end(); iA++)
        {
            if ((int)iA->coords.size() !=0)
            {
                iA->coords.clear();
            }
            
            iA->coords.push_back(2.0*(double)rand()/RAND_MAX);
            iA->coords.push_back(2.0*(double)rand()/RAND_MAX);
            iA->coords.push_back(2.0*(double)rand()/RAND_MAX);        
        }
        /*
        for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA !=allAtoms.end(); iA++)
        {
            std::cout << "X " << std::setprecision(3) << iA->coords[0] << std::endl
                      << "Y "   << std::setprecision(3) << iA->coords[1] << std::endl
                      << "Z "   << std::setprecision(3) << iA->coords[2] << std::endl;
        }
        */
        
        int nCout  = 100;
        int iCout  = 0;
        int nStep  = 100;
        
        int iStep  = 0;
        
        REAL tLamD = 1.0;
        REAL deltaD = 0.1/nCout;
        
        bool lMoved = true;
        
        int iRound = 0;
        while (iCout < nCout && lMoved)
        {
            while (iStep < nStep && lMoved)
            {
                //std::cout << "iStep " << iStep << std::endl
                //          << " nStep " << nStep << std::endl;
                REAL mu = (double) rand()/RAND_MAX;
                if (mu < 0.9)
                {
                    iRound++;
                    // std::cout << "The round " << iRound << std::endl << std::endl;
                    
                    lMoved = false;
                    for (int i=0; i < (int)allAtoms.size(); i++)
                    {
                     
                        for (int j=i+1; j < (int)allAtoms.size(); j++)
                        {
                            //distances by bonds
                            
                            if (fabs(tDistMat[i][j]) > 0.001)
                            {
                                // tDistMat[i][j] > 0.001 is just further protection 
                                REAL tDist, tDiff;
                                tDist = distanceV(allAtoms[i].coords, allAtoms[j].coords);
                                if (tDist < tDistMat[j][i])
                                {
                                    // lower bound
                                    tDiff = (tDistMat[j][i]-tDist)/(tDist + 1.0e-8);
                                }
                                else if (tDist >  tDistMat[i][j])
                                {
                                    // upper bound 
                                    tDiff = (tDistMat[i][j]-tDist)/(tDist + 1.0e-8);
                                }
                                else
                                {
                                    REAL tDiff1 = (tDistMat[i][j] - tDistMat[j][i])/2.0;
                                    tDiff = tDiff1/tDist;
                                }
                                //std::cout << "Diff " << fabs(tDiff) << std::endl;
                                if (fabs(tDiff) > 0.00001)
                                {
                                    for (int k=0; k < (int)allAtoms[i].coords.size(); k++)
                                    {
                                        // std::cout << "tLamD " << tLamD << std::endl;
                                        // std::cout << "increased " << 0.5*tLamD*tDiff*(allAtoms[i].coords[k]-allAtoms[j].coords[k])
                                        //          << std::endl;
                                        allAtoms[i].coords[k]  +=(0.5*tLamD*tDiff*(allAtoms[i].coords[k]-allAtoms[j].coords[k]));
                                        allAtoms[j].coords[k]  +=(0.5*tLamD*tDiff*(allAtoms[j].coords[k]-allAtoms[i].coords[k]));
                                        //std::cout << " atom " << i << " coordinate " << j 
                                        //          << " value " <<  allAtoms[*iNB].coords[j] << std::endl;
                                        // std::cout << " atom " << *iNB << " coordinate " << j 
                                        //          << " value " <<  allAtoms[*iNB].coords[j] << std::endl;
                                        
                                    }
                                    tDist = distanceV(allAtoms[i].coords, allAtoms[j].coords);
                                    tDiff = (tDistMat[i][j]-tDist)/(tDist + 1.0e-8);
                                    //std::cout << "updated Diff " << fabs(tDiff) << std::endl;
                                    lMoved = true;
                                }
                            }
                        }
                    }
                    
                }
                else
                {
                    // should do check_chiral here
                }
                iStep++; 
            }
            tLamD -=deltaD;
            iCout++;
        }
        
        // Check 
        //std::cout << "enter here " << std::endl;
        
        for (int i=0; i <(int)allAtoms.size(); i++)
        {   
            for (std::vector<int>::iterator iNB=allAtoms[i].connAtoms.begin();
                                iNB !=allAtoms[i].connAtoms.end(); iNB++)
            {
                //std::cout << "Atom coordinates for atom " << allAtoms[i].id << std::endl
                //          << "X: " << allAtoms[i].coords[0] << "  Y:  " 
                //          << allAtoms[i].coords[1] << "  Z: " 
                //          << allAtoms[i].coords[2] << std::endl;
                //std::cout << "Atom coordinates for NB atom " << allAtoms[*iNB].id << std::endl
                //          << "X: " << allAtoms[*iNB].coords[0] << "  Y:  " 
                //          << allAtoms[*iNB].coords[1] << "  Z: " 
                //          << allAtoms[*iNB].coords[2] << std::endl;
                REAL tDist = distanceV(allAtoms[i].coords, allAtoms[*iNB].coords);
                std::cout << "Bond length "  << tDistMat[i][*iNB] 
                          << " :  initial distance " << tDist << std::endl;
        
            }
        
        }
        
        for (std::vector<AngleDict>::iterator iAN = allAngles.begin();
                iAN != allAngles.end(); iAN++)
        {
            std::vector<REAL> v21, v31;
            for (int j=0; j < (int)allAtoms[iAN->atoms[0]].coords.size(); j++)
            {
                v21.push_back((allAtoms[iAN->atoms[1]].coords[j]-allAtoms[iAN->atoms[0]].coords[j]));
                v31.push_back((allAtoms[iAN->atoms[2]].coords[j]-allAtoms[iAN->atoms[0]].coords[j]));
            }
            REAL tAngle = getAngle2V(v21, v31);
            std::cout << "The angle between " << allAtoms[iAN->atoms[1]].id << " "
                               << allAtoms[iAN->atoms[0]].id << "  " << allAtoms[iAN->atoms[2]].id 
                               << " : "  << std::endl
                               << " Dictionary value " << iAN->value 
                               << " : initial value "    << tAngle*PID180 << std::endl;
        }
                              
    }
    
    void CodClassify::initRoughCoordsTor()
    {
       
    }
    void CodClassify::outRestraintCif(FileName tFileName)
    {  
        
        std::ofstream outRestrF(tFileName);
        std::string ttNameStr(tFileName);
        std::string tMonoName = ttNameStr.substr(0, 3);
        StrUpper(tMonoName);
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
                    << tMonoName <<"\t"<< tMonoName << "\t" << "'.\t\t'\t"
                    << "non-polymer\t" << (int)allAtoms.size() << "\t.\t."
                    << std::endl;
            
            outRestrF <<"# ------------------------------------------------------" << std::endl
                      <<"# ------------------------------------------------------" << std::endl
                      <<"#" << std::endl
                      <<"# --- DESCRIPTION OF MONOMERS ---" << std::endl
                      <<"#" << std::endl
                      <<"data_comp_" << tMonoName << std::endl
                      <<"#" << std::endl; 
                    
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
            
            for (std::vector<AtomDict>::iterator iA = allAtoms.begin();
                    iA != allAtoms.end(); iA++)
            {
                //double r1 =  (double) rand()/RAND_MAX;
                //double r2 =  (double) rand()/RAND_MAX;
                //double r3 =  (double) rand()/RAND_MAX;
                outRestrF << tMonoName << "\t" 
                          << iA->id << "\t"
                          << iA->chemType << "\t"
                          << iA->enerType << "\t" 
                          << iA->parCharge << "\t"
                          << std::setprecision(3) << std::fixed 
                          << iA->coords[0] << "\t" 
                          << std::setprecision(3) << std::fixed 
                          << iA->coords[1] << "\t"
                          << std::setprecision(3) << std::fixed 
                          << iA->coords[2] << std::endl;
            }
            
            // Bond sections 
            outRestrF << "loop_" << std::endl
                    << "_chem_comp_bond.comp_id" << std::endl
                    << "_chem_comp_bond.atom_id_1" << std::endl
                    << "_chem_comp_bond.atom_id_2" << std::endl
                    << "_chem_comp_bond.type" << std::endl
                    << "_chem_comp_bond.value_dist"<< std::endl
                    << "_chem_comp_bond.value_dist_esd" << std::endl;
                   // << "_chem_comp_bond.exact_cod_dist" << std::endl;
            
            for (std::vector<BondDict>::iterator iB=allBonds.begin();
                    iB !=allBonds.end(); iB++)
            {
                outRestrF << tMonoName << "\t" 
                          << iB->atoms[0] << "\t"
                          << iB->atoms[1] << "\t"
                          << iB->order <<"\t"
                          << std::setprecision(3) << std::fixed
                          << iB->value << "\t" 
                          << std::setprecision(3) << std::fixed
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
            
            // Angle section
            outRestrF << "loop_" << std::endl
                      << "_chem_comp_angle.comp_id"   << std::endl
                      << "_chem_comp_angle.atom_id_1" << std::endl
                      << "_chem_comp_angle.atom_id_2" << std::endl
                      << "_chem_comp_angle.atom_id_3" << std::endl
                      << "_chem_comp_angle.value_angle"     << std::endl
                      << "_chem_comp_angle.value_angle_esd" << std::endl;
                    //  << "_chem_comp_angle.exact_cod_dist"  << std::endl;
            
            for (std::vector<AngleDict>::iterator iA=allAngles.begin();
                    iA != allAngles.end(); iA++)
            {
                //for (std::vector<int>::iterator iAt=iA->atoms.begin();
                //        iAt !=iA->atoms.end(); iAt++)
                //{
                // difference in comp atom definitions between cod and
                // dictionary: inner-out1-out2(cod),
                // atom1-atom2(center)-atom3(dictionary)
                
                outRestrF << tMonoName << "\t"
                          << allAtoms[iA->atoms[1]].id << "\t"
                          << allAtoms[iA->atoms[0]].id << "\t"
                          << allAtoms[iA->atoms[2]].id << "\t";
                
                //}
                        
                outRestrF << std::setprecision(3) << iA->value << "\t"
                        << std::setprecision(2) << iA->sigValue << std::endl;
                
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
            
            // Torsion section 
            
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
            for (std::vector<TorsionDict>::iterator iT=allTorsions.begin();
                    iT !=allTorsions.end(); iT++)
            {
                // std::string idxTorStr=IntToStr(idxTor);
                // idxTorStr = "tor_" + idxTorStr;
                outRestrF << tMonoName << "\t" 
                        << iT->id <<"\t" 
                        <<  allAtoms[iT->atoms[0]].id <<"\t"
                        <<  allAtoms[iT->atoms[1]].id <<"\t"
                        <<  allAtoms[iT->atoms[2]].id <<"\t"
                        <<  allAtoms[iT->atoms[3]].id <<"\t"
                        << iT->value << "\t" 
                        << std::setprecision(2) << std::fixed << "10.00\t" 
                        << iT->period << std::endl;
                idxTor++;        
            }
            
            //  For chiral centers
            outRestrF << "loop_" << std::endl
                    << "_chem_comp_chir.comp_id" << std::endl
                    << "_chem_comp_chir.id" << std::endl
                    << "_chem_comp_chir.atom_id_centre" << std::endl
                    << "_chem_comp_chir.atom_id_1" << std::endl
                    << "_chem_comp_chir.atom_id_2" << std::endl
                    << "_chem_comp_chir.atom_id_3" << std::endl
                    << "_chem_comp_chir.volume_sign" << std::endl;
            
            int idxC =1;
            for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                    iA !=allAtoms.end(); iA++)
            {
                if (iA->chiralIdx == 1)
                {
                    std::vector<ID> chirAtms;
                    for (std::vector<int>::iterator iNB=iA->connAtoms.begin();
                            iNB != iA->connAtoms.end(); iNB++)
                    {
                        std::size_t tFind = allAtoms[*iNB].chemType.find("H");
                        if (tFind==std::string::npos)
                        {
                            chirAtms.push_back(allAtoms[*iNB].id);
                        }
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
                        outRestrF << tMonoName << "\t" 
                                  << idxStr << "\t"
                                  << iA->id  << "\t" ;
                       
                        outRestrF << chirAtms[0] << "\t";
                        outRestrF << chirAtms[1] << "\t";
                        outRestrF << chirAtms[2] << "\t";
                        outRestrF << "BOTH" << std::endl;
                        idxC++;
                    }
                }
            }
            
            
            // Planar group section
            outRestrF << "loop_" << std::endl
                    << "_chem_comp_plane_atom.comp_id"  << std::endl
                    << "_chem_comp_plane_atom.plane_id" << std::endl
                    << "_chem_comp_plane_atom.atom_id"  << std::endl
                    << "_chem_comp_plane_atom.dist_esd" << std::endl;
            
            int idxP = 1;
            for (std::vector<PlaneDict>::iterator iP=allPlanes.begin();
                    iP !=allPlanes.end(); iP++)
            {
                std::string idxPStr = IntToStr(idxP);
                idxPStr = "plan-" + idxPStr;
                for(std::map<ID, int>::iterator iAt=iP->atoms.begin();
                        iAt != iP->atoms.end(); iAt++)
                {
                    outRestrF << tMonoName << "\t" << idxPStr + "\t"
                            << iAt->first << "\t" 
                            << "0.020\t" << std::endl;
                }
                idxP++;
            }
                                 
            outRestrF.close();
            
        }
    }
    
    void CodClassify::outRestraintCif(FileName tFileName, ID tMonoName)
    {  
        std::ofstream outRestrF(tFileName);
        
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
                    << tMonoName <<"\t"<< tMonoName << "\t" << "'.\t\t'\t"
                    << "non-polymer\t" << (int)allAtoms.size() << "\t" 
                    << (int)allHydroAtoms.size() << "\t."
                    << std::endl;
            
            outRestrF <<"# ------------------------------------------------------" << std::endl
                      <<"# ------------------------------------------------------" << std::endl
                      <<"#" << std::endl
                      <<"# --- DESCRIPTION OF MONOMERS ---" << std::endl
                      <<"#" << std::endl
                      <<"data_comp_" << tMonoName << std::endl
                      <<"#" << std::endl; 
                    
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
            
            for (std::vector<AtomDict>::iterator iA = allAtoms.begin();
                    iA != allAtoms.end(); iA++)
            {
                //double r1 =  (double) rand()/RAND_MAX;
                //double r2 =  (double) rand()/RAND_MAX;
                //double r3 =  (double) rand()/RAND_MAX;
                StrUpper(iA->chemType);
                outRestrF << tMonoName << "\t" 
                          << iA->id << "\t"
                          << iA->chemType << "\t"
                          << iA->enerType << "\t" 
                          << iA->parCharge << "\t"
                          << std::setprecision(3) << std::fixed 
                          << iA->coords[0] << "\t" 
                          << std::setprecision(3) << std::fixed 
                          << iA->coords[1] << "\t"
                          << std::setprecision(3) << std::fixed 
                          << iA->coords[2] << std::endl;
            }
            
            
            // Bond sections 
            outRestrF << "loop_" << std::endl
                    << "_chem_comp_bond.comp_id" << std::endl
                    << "_chem_comp_bond.atom_id_1" << std::endl
                    << "_chem_comp_bond.atom_id_2" << std::endl
                    << "_chem_comp_bond.type" << std::endl
                    << "_chem_comp_bond.value_dist"<< std::endl
                    << "_chem_comp_bond.value_dist_esd" << std::endl;
                   // << "_chem_comp_bond.exact_cod_dist" << std::endl;
            
            for (std::vector<BondDict>::iterator iB=allBonds.begin();
                    iB !=allBonds.end(); iB++)
            {
                outRestrF << tMonoName << "\t" 
                          << iB->atoms[0] << "\t"
                          << iB->atoms[1] << "\t"
                          << iB->order << "\t"
                          << std::setprecision(3) << std::fixed
                          << iB->value << "\t" 
                          << std::setprecision(3) << std::fixed
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
           
            // Angle section
            outRestrF << "loop_" << std::endl
                      << "_chem_comp_angle.comp_id"   << std::endl
                      << "_chem_comp_angle.atom_id_1" << std::endl
                      << "_chem_comp_angle.atom_id_2" << std::endl
                      << "_chem_comp_angle.atom_id_3" << std::endl
                      << "_chem_comp_angle.value_angle"     << std::endl
                      << "_chem_comp_angle.value_angle_esd" << std::endl;
                    //  << "_chem_comp_angle.exact_cod_dist"  << std::endl;
            
            for (std::vector<AngleDict>::iterator iA=allAngles.begin();
                    iA != allAngles.end(); iA++)
            {
                //for (std::vector<int>::iterator iAt=iA->atoms.begin();
                //        iAt !=iA->atoms.end(); iAt++)
                //{
                // difference in comp atom definitions between cod and
                // dictionary: inner-out1-out2(cod),
                // atom1-atom2(center)-atom3(dictionary)
              
                if (allAtoms[iA->atoms[0]].isMetal)
                {
                    for (std::vector<REAL>::iterator iMA=iA->codAngleValues.begin();
                            iMA != iA->codAngleValues.end(); iMA++)
                    {
                            outRestrF << tMonoName << "\t"
                                      << allAtoms[iA->atoms[1]].id << "\t"
                                      << allAtoms[iA->atoms[0]].id << "\t"
                                      << allAtoms[iA->atoms[2]].id << "\t";
                            outRestrF << std::setprecision(3) << *iMA << "\t"
                                      <<  3.00 << std::endl;
                    }
                }
                else
                {
                        outRestrF << tMonoName << "\t"
                              << allAtoms[iA->atoms[1]].id << "\t"
                              << allAtoms[iA->atoms[0]].id << "\t"
                              << allAtoms[iA->atoms[2]].id << "\t";
                        outRestrF << std::setprecision(3) << iA->value << "\t"
                             << std::setprecision(2) << iA->sigValue << std::endl;
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
            
            // Torsion section 
            if((int)allTorsions.size() !=0)
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
            
                for (std::vector<TorsionDict>::iterator iT=allTorsions.begin();
                        iT !=allTorsions.end(); iT++)
                {
                    //std::string idxTorStr=IntToStr(idxTor);
                    //idxTorStr = "tor_" + idxTorStr;
                    // std::cout << "Torsion angle " << idxTor 
                    //          << " It contains " << (int)iT->atoms.size() << std::endl;
                          
                    //std::cout << iT->atoms[0] << std::endl
                    //          << iT->atoms[1] << std::endl
                    //          << iT->atoms[2] << std::endl
                    //          << iT->atoms[3] << std::endl;
                
                    outRestrF << tMonoName << "\t" 
                            << iT->id  <<"\t" 
                            <<  allAtoms[iT->atoms[0]].id <<"\t"
                            <<  allAtoms[iT->atoms[1]].id <<"\t"
                            <<  allAtoms[iT->atoms[2]].id <<"\t"
                            <<  allAtoms[iT->atoms[3]].id <<"\t"
                            << iT->value << "\t" 
                            << std::setprecision(2) << std::fixed << "10.00\t" 
                            << iT->period << std::endl;
                    idxTor++;        
                }
                
            }
            
            //  For chiral centers
            bool l_ch = false;
            if ((int)allChirals.size() !=0)
            {
                l_ch = true;
            }
            else
            {
                for (int i_ch =0; i_ch < (int)allAtoms.size(); i_ch++)
                {
                    if (allAtoms[i_ch].chiralIdx ==1)
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
                for (std::vector<ChiralDict>::iterator iCh = allChirals.begin();
                        iCh != allChirals.end(); iCh++)
                {
                    inputChiralID.push_back(iCh->archID);
                    outRestrF << tMonoName << "\t" 
                              << iCh->id  << "\t";
                    for (std::vector<int>::iterator iAt=iCh->atoms.begin();
                            iAt != iCh->atoms.end(); iAt++)
                    {
                        outRestrF << allAtoms[*iAt].id << "\t";
                    }
                    outRestrF << iCh->sign << std::endl;
                }
                // New chiral that are not in the input list 
                int idxC =(int)inputChiralID.size();
                for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                          iA !=allAtoms.end(); iA++)
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
                                
                                if (allAtoms[*iNB].chemType !="H")
                                {
                                    chirAtms.push_back(allAtoms[*iNB].id);
                                }
                                else
                                {
                                    HAtms.push_back(allAtoms[*iNB].id);
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
                                  
                                outRestrF << tMonoName << "\t" 
                                      << idxStr << "\t"
                                      << iA->id  << "\t" ;
                       
                                outRestrF << chirAtms[0] << "\t";
                                outRestrF << chirAtms[1] << "\t";
                                outRestrF << chirAtms[2] << "\t";
                                outRestrF << "BOTH" << std::endl;
                                idxC++;
                            }
                        }
                    }
                }
            }
            
            
            
            // Planar group section
            if ((int)allPlanes.size() >0)
            {
                outRestrF << "loop_" << std::endl
                          << "_chem_comp_plane_atom.comp_id"  << std::endl
                          << "_chem_comp_plane_atom.plane_id" << std::endl
                          << "_chem_comp_plane_atom.atom_id"  << std::endl
                          << "_chem_comp_plane_atom.dist_esd" << std::endl;
                int idxP = 1;
                for (std::vector<PlaneDict>::iterator iP=allPlanes.begin();
                        iP !=allPlanes.end(); iP++)
                {
                    if ((int)iP->atoms.size() >3)
                    {
                        std::string idxPStr = IntToStr(idxP);
                        idxPStr = "plan-" + idxPStr;
                        for(std::map<ID, int>::iterator iAt=iP->atoms.begin();
                               iAt != iP->atoms.end(); iAt++)
                        {
                            outRestrF << tMonoName << "\t" << idxPStr + "\t"
                                      << iAt->first << "\t" 
                                      << "0.020\t" << std::endl;
                        }
                        idxP++;
                    }
                }
            }
            outRestrF.close();
        }
    }
    
 
    void CodClassify::outRestraintCif2(FileName tFileName, ID tMonoName)
    {  
        std::ofstream outRestrF(tFileName);
        
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
                    << tMonoName <<"\t"<< tMonoName << "\t" << "'.\t\t'\t"
                    << "non-polymer\t" << (int)allAtoms.size() << "\t" 
                    << (int)allHydroAtoms.size() << "\t."
                    << std::endl;
            
            outRestrF <<"# ------------------------------------------------------" << std::endl
                      <<"# ------------------------------------------------------" << std::endl
                      <<"#" << std::endl
                      <<"# --- DESCRIPTION OF MONOMERS ---" << std::endl
                      <<"#" << std::endl
                      <<"data_comp_" << tMonoName << std::endl
                      <<"#" << std::endl; 
                    
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
            
            for (std::vector<AtomDict>::iterator iA = allAtoms.begin();
                    iA != allAtoms.end(); iA++)
            {
                //double r1 =  (double) rand()/RAND_MAX;
                //double r2 =  (double) rand()/RAND_MAX;
                //double r3 =  (double) rand()/RAND_MAX;
                StrUpper(iA->chemType);
                outRestrF << tMonoName << "\t" 
                          << iA->id << "\t"
                          << iA->chemType << "\t"
                          << iA->enerType << "\t" 
                          << iA->parCharge << "\t"
                          << std::setprecision(3) << std::fixed 
                          << iA->coords[0] << "\t" 
                          << std::setprecision(3) << std::fixed 
                          << iA->coords[1] << "\t"
                          << std::setprecision(3) << std::fixed 
                          << iA->coords[2] << std::endl;
            }
            
            
            // Bond sections 
            outRestrF << "loop_" << std::endl
                    << "_chem_comp_bond.comp_id" << std::endl
                    << "_chem_comp_bond.atom_id_1" << std::endl
                    << "_chem_comp_bond.atom_id_2" << std::endl
                    << "_chem_comp_bond.type" << std::endl
                    << "_chem_comp_bond.value_dist"<< std::endl
                    << "_chem_comp_bond.value_dist_esd" << std::endl;
                   // << "_chem_comp_bond.exact_cod_dist" << std::endl;
            
            for (std::vector<BondDict>::iterator iB=allBonds.begin();
                    iB !=allBonds.end(); iB++)
            {
                outRestrF << tMonoName << "\t" 
                          << iB->atoms[0] << "\t"
                          << iB->atoms[1] << "\t"
                          << iB->order << "\t"
                          << std::setprecision(3) << std::fixed
                          << iB->value << "\t" 
                          << std::setprecision(3) << std::fixed
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
           
            // Angle section
            outRestrF << "loop_" << std::endl
                      << "_chem_comp_angle.comp_id"   << std::endl
                      << "_chem_comp_angle.atom_id_1" << std::endl
                      << "_chem_comp_angle.atom_id_2" << std::endl
                      << "_chem_comp_angle.atom_id_3" << std::endl
                      << "_chem_comp_angle.value_angle"     << std::endl
                      << "_chem_comp_angle.value_angle_esd" << std::endl;
                    //  << "_chem_comp_angle.exact_cod_dist"  << std::endl;
            
            for (std::vector<AngleDict>::iterator iA=allAngles.begin();
                    iA != allAngles.end(); iA++)
            {
                //for (std::vector<int>::iterator iAt=iA->atoms.begin();
                //        iAt !=iA->atoms.end(); iAt++)
                //{
                // difference in comp atom definitions between cod and
                // dictionary: inner-out1-out2(cod),
                // atom1-atom2(center)-atom3(dictionary)
              
                
                outRestrF << tMonoName << "\t"
                          << allAtoms[iA->atoms[1]].id << "\t"
                          << allAtoms[iA->atoms[0]].id << "\t"
                          << allAtoms[iA->atoms[2]].id << "\t";
                outRestrF << std::setprecision(3) << iA->value << "\t"
                          << std::setprecision(2) << iA->sigValue << std::endl;
            }
              
            
            // Torsion section 
            if((int)allTorsions.size() !=0)
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
            
                for (std::vector<TorsionDict>::iterator iT=allTorsions.begin();
                        iT !=allTorsions.end(); iT++)
                {
                    std::string idxTorStr=IntToStr(idxTor);
                    idxTorStr = "tor_" + idxTorStr;
                    // std::cout << "Torsion angle " << idxTor 
                    //          << " It contains " << (int)iT->atoms.size() << std::endl;
                          
                    //std::cout << iT->atoms[0] << std::endl
                    //          << iT->atoms[1] << std::endl
                    //          << iT->atoms[2] << std::endl
                    //          << iT->atoms[3] << std::endl;
                
                    outRestrF << tMonoName << "\t" 
                            << idxTorStr <<"\t" 
                            <<  allAtoms[iT->atoms[0]].id <<"\t"
                            <<  allAtoms[iT->atoms[1]].id <<"\t"
                            <<  allAtoms[iT->atoms[2]].id <<"\t"
                            <<  allAtoms[iT->atoms[3]].id <<"\t"
                            << iT->value << "\t" 
                            << std::setprecision(2) << std::fixed << "10.00\t" 
                            << iT->period << std::endl;
                    idxTor++;        
                }
                
            }
            
            //  For chiral centers
            bool l_ch = false;
            if ((int)allChirals.size() !=0)
            {
                l_ch = true;
            }
            else
            {
                for (int i_ch =0; i_ch < (int)allAtoms.size(); i_ch++)
                {
                    if (allAtoms[i_ch].chiralIdx ==1)
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
                for (std::vector<ChiralDict>::iterator iCh = allChirals.begin();
                        iCh != allChirals.end(); iCh++)
                {
                    inputChiralID.push_back(iCh->archID);
                    outRestrF << tMonoName << "\t" 
                              << iCh->id  << "\t";
                    for (std::vector<int>::iterator iAt=iCh->atoms.begin();
                            iAt != iCh->atoms.end(); iAt++)
                    {
                        outRestrF << allAtoms[*iAt].id << "\t";
                    }
                    outRestrF << iCh->sign << std::endl;
                }
                // New chiral that are not in the input list 
                int idxC =(int)inputChiralID.size();
                for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                          iA !=allAtoms.end(); iA++)
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
                                
                                if (allAtoms[*iNB].chemType !="H")
                                {
                                    chirAtms.push_back(allAtoms[*iNB].id);
                                }
                                else
                                {
                                    HAtms.push_back(allAtoms[*iNB].id);
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
                                  
                                outRestrF << tMonoName << "\t" 
                                      << idxStr << "\t"
                                      << iA->id  << "\t" ;
                       
                                outRestrF << chirAtms[0] << "\t";
                                outRestrF << chirAtms[1] << "\t";
                                outRestrF << chirAtms[2] << "\t";
                                outRestrF << "BOTH" << std::endl;
                                idxC++;
                            }
                        }
                    }
                }
            }
            
            
            
            // Planar group section
            if ((int)allPlanes.size() >0)
            {
                outRestrF << "loop_" << std::endl
                          << "_chem_comp_plane_atom.comp_id"  << std::endl
                          << "_chem_comp_plane_atom.plane_id" << std::endl
                          << "_chem_comp_plane_atom.atom_id"  << std::endl
                          << "_chem_comp_plane_atom.dist_esd" << std::endl;
                int idxP = 1;
                for (std::vector<PlaneDict>::iterator iP=allPlanes.begin();
                        iP !=allPlanes.end(); iP++)
                {
                    if ((int)iP->atoms.size() >3 )
                    {
                        std::string idxPStr = IntToStr(idxP);
                        idxPStr = "plan-" + idxPStr;
                        for(std::map<ID, int>::iterator iAt=iP->atoms.begin();
                                   iAt != iP->atoms.end(); iAt++)
                        {
                            outRestrF << tMonoName << "\t" << idxPStr + "\t"
                                      << iAt->first << "\t" 
                                     << "0.020\t" << std::endl;
                        }
                        idxP++;
                    }
                }
            }
            outRestrF.close();
        }
    }
    
    void CodClassify::outPDB(FileName tFName)
    {
        // This is a temporary one, the method should be defined outside 
        // this class.
        std::string tName(tFName);
        std::string ttNameStr(tFName);
        std::string tMonoName = ttNameStr.substr(0, 3);
        StrUpper(tMonoName);
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
            outPDB << std::left << "MONOMER tests " << tMonoName;
            outPDB.width(10);
            outPDB << std::left << "Date:";
            outPDB.width(6);
            outPDB << std::left << tMonoName;
            outPDB.width(14);
            outPDB << std::left <<"" <<std::endl;
                
            // CRYST1 section 
            outPDB.width(80);
            outPDB << std::left 
                   <<"CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1" 
                   <<std::endl;
            
            // ATOM sections
            for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                    iA !=allAtoms.end(); iA++)
            {
                //double r1 =  (double) rand()/RAND_MAX;
                //double r2 =  (double) rand()/RAND_MAX;
                //double r3 =  (double) rand()/RAND_MAX;
                
                outPDB.width(6);
                outPDB <<std::left<< "ATOM";
                outPDB.width(5);
                outPDB <<std::right << iA->seriNum;
                outPDB.width(1);
                outPDB << " ";
                if ((int)iA->id.size() ==4)
                {
                    outPDB.width(4);
                    outPDB << std::left << iA->id;
                }
                else
                {
                    outPDB.width(1);
                    outPDB << " ";
                    outPDB.width(3);
                    outPDB << std::left << iA->id;
                }
                outPDB.width(1); // altLoc
                outPDB << " ";
                outPDB.width(3); // resName
                outPDB << "XXX";
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

    void CodClassify::outPDB(FileName tFName, ID tMonoName)
    {
        // This is a temporary one, the method should be defined outside 
        // this class.
        std::string tName(tFName);    
        StrUpper(tMonoName);
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
            outPDB << std::left << " MONOMER tests " << tMonoName;
            outPDB.width(10);
            outPDB << std::left << " Date:";
            outPDB.width(6);
            outPDB << std::left << tMonoName;
            outPDB.width(14);
            outPDB << std::left <<"" <<std::endl;
                
            // CRYST1 section 
            outPDB.width(80);
            outPDB << std::left 
                   <<"CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1" 
                   <<std::endl;
            
            // ATOM sections
            for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                    iA !=allAtoms.end(); iA++)
            {
                //double r1 =  (double) rand()/RAND_MAX;
                //double r2 =  (double) rand()/RAND_MAX;
                //double r3 =  (double) rand()/RAND_MAX;
                
                outPDB.width(6);
                outPDB <<std::left<< "ATOM";
                outPDB.width(5);
                outPDB <<std::right << iA->seriNum;
                outPDB.width(1);
                outPDB << " ";
                if ((int)iA->id.size() ==4)
                {
                    outPDB.width(4);
                    outPDB << std::left << iA->id;
                }
                else
                {
                    outPDB.width(1);
                    outPDB << " ";
                    outPDB.width(3);
                    outPDB << std::left << iA->id;
                }
                outPDB.width(1); // altLoc
                outPDB << " ";
                outPDB.width(3); // resName
                outPDB << tMonoName.substr(0,3);
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
    
    void CodClassify::outAtomTypes(ID tMonoRootName)
    {
        std::string    outAtomTypeName = tMonoRootName + "_CodAtomType.txt";
        std::ofstream  outAtomtype(outAtomTypeName.c_str());
        
        if (outAtomtype.is_open())
        {
            outAtomtype << "Number of Atoms in the file " 
                        << (int)allAtoms.size() << std::endl;
            
            for (std::vector<AtomDict>::iterator iA=allAtoms.begin();
                    iA !=allAtoms.end(); iA++)
            {
                outAtomtype << iA->seriNum+1  << "\t"
                            << iA->id << "\t" 
                            << iA->codClass << std::endl; 
            }
            
            outAtomtype.close(); 
        }
        
    }
    
    void CodClassify::getAnglesFromPDB(ID tFName)
    {
        DictPDBFile refPDB(tFName, std::ios::in);
        
        // put atomic coordinates after refinement into atoms of AtomDict
        for(std::vector<AtomDict>::iterator iA=allAtoms.begin();
                iA != allAtoms.end(); iA++)
        {
            std::map<ID, std::vector<REAL> >::iterator iFind;
            iFind = refPDB.allHetAtmList.find(iA->id);
            if (iFind != refPDB.allHetAtmList.end())
            {
                iA->coords[0] = refPDB.allHetAtmList[iA->id][0];
                iA->coords[1] = refPDB.allHetAtmList[iA->id][1];
                iA->coords[2] = refPDB.allHetAtmList[iA->id][2];
            }
            else
            {
                std::cout << "Could not find coordinate values for atom " 
                        << iA->id << std::endl;
                exit(1);
            }
        }
        
        // Now calculate unique angles among atoms, some of which may have 
        // several candidate values
        std::map<ID,int> angleIDs;
        for (std::vector<AngleDict>::iterator iAn=allAngles.begin();
                    iAn != allAngles.end(); iAn++)
        {
            std::vector<REAL> aV1, aV2;
            for (int i=0; i < (int)allAtoms[iAn->atoms[0]].coords.size(); i++)
            {
                aV1.push_back(allAtoms[iAn->atoms[1]].coords[i]-allAtoms[iAn->atoms[0]].coords[i]);
                aV2.push_back(allAtoms[iAn->atoms[2]].coords[i]-allAtoms[iAn->atoms[0]].coords[i]);
            }
            iAn->value = getAngle2V(aV1, aV2)*PID180;
        }
        
    }
    
    void CodClassify::setSpecial3NBSymb(std::vector<AtomDict>::iterator tAt)
    {
        
        std::vector<int> serNumNB12;
        std::map<std::string, int>   NB3Props;
        
        if (tAt->ringRep.size() !=0)
        {
            for (std::vector<int>::iterator iNB1=tAt->connAtoms.begin();
                  iNB1 !=tAt->connAtoms.end(); iNB1++)
            {
                if (std::find(serNumNB12.begin(), serNumNB12.end(), *iNB1)
                        == serNumNB12.end())
                {
                    serNumNB12.push_back(*iNB1);
                }
                for (std::vector<int>::iterator iNB2=allAtoms[*iNB1].connAtoms.begin();
                            iNB2!=allAtoms[*iNB1].connAtoms.end(); iNB2++)
                {
                    if (std::find(serNumNB12.begin(), serNumNB12.end(), *iNB2)
                        == serNumNB12.end())
                    {
                        serNumNB12.push_back(*iNB2);
                    }
                }
            }
                
            
            for (std::vector<int>::iterator iNB1=tAt->connAtoms.begin();
                    iNB1 !=tAt->connAtoms.end(); iNB1++)
            {
                if (allAtoms[*iNB1].ringRep.size() !=0)
                {
                    for (std::vector<int>::iterator iNB2=allAtoms[*iNB1].connAtoms.begin();
                          iNB2 !=allAtoms[*iNB1].connAtoms.end(); iNB2++)
                    {
                        if (allAtoms[*iNB2].ringRep.size() !=0)
                        {
                            for (std::vector<int>::iterator iNB3=allAtoms[*iNB2].connAtoms.begin();
                                 iNB3 !=allAtoms[*iNB2].connAtoms.end(); iNB3++)
                            {
                                if (std::find(serNumNB12.begin(), serNumNB12.end(), *iNB3)
                                      ==serNumNB12.end())
                                {
                                    
                                    std::string tProp = allAtoms[*iNB3].chemType; 
                                    tProp.append("<");        
                                    tProp.append(IntToStr((int)allAtoms[*iNB3].connAtoms.size()));
                                    tProp.append(">");                    
                                                       
                                    if (NB3Props.find(tProp)==NB3Props.end())
                                    {
                                        NB3Props[tProp] = 1;
                                    }
                                    else
                                    {
                                        NB3Props[tProp]++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            
            // Now symbols 
            std::list<std::string> tComps;
            
            for (std::map<std::string, int>::iterator iNB3=NB3Props.begin();
                    iNB3 !=NB3Props.end(); iNB3++)
            {
                std::string tID = "";
                tID.append(IntToStr(iNB3->second));
                tID.append("|");
                tID.append(iNB3->first);
                tComps.push_back(tID);
            }
            
            tComps.sort(compareNoCase2);
            
            std::cout << "The following are special 3 NB around atom " 
                      << tAt->id << std::endl;
            
            for (std::list<std::string>::iterator iID=tComps.begin();
                    iID !=tComps.end(); iID++)
            {
                std::cout << *iID << std::endl;
            }
            
            
            
            ID all3 = "{";
            
            unsigned i=0, aN=tComps.size();
            
            for (std::list<std::string>::iterator iID=tComps.begin();
                    iID !=tComps.end(); iID++)
            {
                if (i < aN-1)
                {
                    all3.append(*iID + ",");
                }
                else
                {
                    all3.append(*iID);
                }
                i++;
            }
            
            all3.append("}");
        
            std::cout << "The special 3 NB string is " << all3 << std::endl;
            
            tAt->codClass.append(all3);
            
        }
        
    }
    
    void CodClassify::setSpecial3NBSymb2(std::vector<AtomDict>::iterator tAt)
    {
        
        std::cout << "For atom " << tAt->id 
                  << " of serial number " << tAt->seriNum << " : " << std::endl;
        
        std::vector<int> serNumNB123;
        std::map<std::string, int>   NB3Props;
        
        if (tAt->ringRep.size() !=0)
        {
            for (std::vector<int>::iterator iNB1=tAt->connAtoms.begin();
                  iNB1 !=tAt->connAtoms.end(); iNB1++)
            {
                if (std::find(serNumNB123.begin(), serNumNB123.end(), *iNB1)
                        == serNumNB123.end())
                {
                    serNumNB123.push_back(*iNB1);
                }
                for (std::vector<int>::iterator iNB2=allAtoms[*iNB1].connAtoms.begin();
                            iNB2!=allAtoms[*iNB1].connAtoms.end(); iNB2++)
                {
                    if (std::find(serNumNB123.begin(), serNumNB123.end(), *iNB2)
                        == serNumNB123.end())
                    {
                        serNumNB123.push_back(*iNB2);
                    }
                }
            }
                
            
            for (std::vector<int>::iterator iNB1=tAt->connAtoms.begin();
                    iNB1 !=tAt->connAtoms.end(); iNB1++)
            {
                if (allAtoms[*iNB1].ringRep.size() !=0)
                {
                    for (std::vector<int>::iterator iNB2=allAtoms[*iNB1].connAtoms.begin();
                          iNB2 !=allAtoms[*iNB1].connAtoms.end(); iNB2++)
                    {
                        if (allAtoms[*iNB2].ringRep.size() !=0)
                        {
                            for (std::vector<int>::iterator iNB3=allAtoms[*iNB2].connAtoms.begin();
                                 iNB3 !=allAtoms[*iNB2].connAtoms.end(); iNB3++)
                            {
                                if (std::find(serNumNB123.begin(), serNumNB123.end(), *iNB3)
                                      ==serNumNB123.end())
                                {
                                    
                                    std::string tProp = allAtoms[*iNB3].chemType; 
                                    tProp.append("<");        
                                    tProp.append(IntToStr((int)allAtoms[*iNB3].connAtoms.size()));
                                    tProp.append(">");                    
                                                       
                                    if (NB3Props.find(tProp)==NB3Props.end())
                                    {
                                        NB3Props[tProp] = 1;
                                    }
                                    else
                                    {
                                        NB3Props[tProp]++;
                                    }
                                    serNumNB123.push_back(*iNB3);
                                    std::cout << "add 3NB atom " << allAtoms[*iNB3].id << std::endl;
                                }
                            }
                        }
                    }
                }
            }
            
            // Now symbols 
            std::list<std::string> tComps;
            
            for (std::map<std::string, int>::iterator iNB3=NB3Props.begin();
                    iNB3 !=NB3Props.end(); iNB3++)
            {
                std::string tID = "";
                tID.append(IntToStr(iNB3->second));
                tID.append("|");
                tID.append(iNB3->first);
                tComps.push_back(tID);
            }
            
            tComps.sort(compareNoCase2);
            
            std::cout << "The following are special 3 NB around atom " 
                      << tAt->id << std::endl;
            
            for (std::list<std::string>::iterator iID=tComps.begin();
                    iID !=tComps.end(); iID++)
            {
                std::cout << *iID << std::endl;
            }
            
            
            
            ID all3 = "{";
            
            unsigned i=0, aN=tComps.size();
            
            for (std::list<std::string>::iterator iID=tComps.begin();
                    iID !=tComps.end(); iID++)
            {
                if (i < aN-1)
                {
                    all3.append(*iID + ",");
                }
                else
                {
                    all3.append(*iID);
                }
                i++;
            }
            
            all3.append("}");
        
            std::cout << "The special 3NB string is " << all3 << std::endl;
            
            tAt->codClass.append(all3);
            
        }
        
    }
    
    
    /* ###################### Class CodBonds  ##################### */
    
    CodBonds::CodBonds()
    {
    }
    
    CodBonds::CodBonds(const CodBonds & tCBonds)
    {
        for (std::vector<BondDict>::const_iterator iCB=tCBonds.allBonds.begin();
                iCB != tCBonds.allBonds.end(); iCB++)
        {
            allBonds.push_back(*iCB);
        }
        
        for (std::map<ID, std::vector<REAL> >::const_iterator iCM=tCBonds.allBondsMap.begin();
                iCM != tCBonds.allBondsMap.end(); iCM++)
        {
            for (std::vector<REAL>::const_iterator iVA = iCM->second.begin();
                    iVA != iCM->second.end(); iVA++)
            {
                allBondsMap[iCM->first].push_back(*iVA);
            }
        }
    }
    
    CodBonds::CodBonds(FileName                tFname,
                       std::ios::openmode      tOpenMode=std::ios::in)
    {
        if (tOpenMode == std::ios::in)
        {
            
            inFile.open(tFname, tOpenMode);
            
            if (inFile.is_open())
            {
                
                setupSystem();
                
                inFile.close();
                
                 
                /*
                // Check
                std::cout << "all cod bonds: " << std::endl;
                for (std::vector<BondDict>::iterator iCB = allBonds.begin();
                        iCB != allBonds.end(); iCB++)
                {
                    std::cout << "BOND " << iCB->seriNum << std::endl;
                    std::cout << "Atom 1: " << iCB->atoms[0] << std::endl
                            << "its COD_class " << iCB->atomsCodClasses[0]
                            << std::endl
                            << "Atom 2: " << iCB->atoms[1] << std::endl
                            << "its COD_class " << iCB->atomsCodClasses[1]
                            << std::endl
                            << "its Value " << iCB->length << std::endl;
                }
                */
            }
        }
    }
    
    CodBonds::~CodBonds()
    {
        if(inFile.is_open())
        {
           inFile.close(); 
        }
    }
    
    void CodBonds::setupSystem()
    {
        std::string tRecord="";
        
        int iLine = 0;
        
        while(!inFile.eof())
        {
            std::getline(inFile, tRecord);
            tRecord = TrimSpaces(tRecord);
            std::vector<std::string> tBuf;
            StrTokenize(tRecord, tBuf);
            
            if ((int)tBuf.size() ==7)
            {
                
                BondDict  aBond;
                
                aBond.seriNum    = iLine;
                iLine++;
                aBond.value     = StrToReal(TrimSpaces(tBuf[0]));

                aBond.atoms.push_back(TrimSpaces(tBuf[1]));
                aBond.atomsCodClasses.push_back(TrimSpaces(tBuf[2]));
                
                aBond.atoms.push_back(TrimSpaces(tBuf[3]));
                aBond.atomsCodClasses.push_back(TrimSpaces(tBuf[4]));
               
                allBonds.push_back(aBond);
            }
        }
    }
    
    void CodBonds::setupSystem(std::vector<sortLine> & tAllLines)
    {
        if((int)tAllLines.size() > 0)
        {
            int iLine = 0;
            std::string tRecord="";
            for (std::vector<sortLine>::iterator iLi=tAllLines.begin();
                    iLi != tAllLines.end(); iLi++)
            {
                tRecord = TrimSpaces(iLi->line);
                std::vector<std::string> tBuf;
                StrTokenize(tRecord, tBuf);
            
                if ((int)tBuf.size() ==7)
                {
                    BondDict  aBond;
                    
                    aBond.seriNum    = iLine;
                    iLine++;
                    
                    aBond.value     = StrToReal(TrimSpaces(tBuf[0]));

                    aBond.atoms.push_back(TrimSpaces(tBuf[1]));
                    aBond.atomsCodClasses.push_back(TrimSpaces(tBuf[2]));
                
                    aBond.atoms.push_back(TrimSpaces(tBuf[3]));
                    aBond.atomsCodClasses.push_back(TrimSpaces(tBuf[4]));
                               
                    allBonds.push_back(aBond);
                }
            }
        }
    }
    
    void CodBonds::sortCodTable(FileName   tFname, 
                                std::vector<sortLine>  & allSortLines)
    {
        inFile.open(tFname, std::ios::in);
        
        if(inFile.is_open())
        {
            
           std::string tRecord="";
           while(!inFile.eof())
           {
               std::getline(inFile, tRecord);
               tRecord = TrimSpaces(tRecord);
               std::vector<std::string> tBuf;
               StrTokenize(tRecord, tBuf);
               
               if ((int)tBuf.size() ==8)
               {
                   sortLine aSL;
                   aSL.key  = tBuf[2].substr(1,(int)tBuf[2].size()-2);
                   aSL.line = tRecord;
                   
                   allSortLines.push_back(aSL);
                   
               }
           }
           
           if((int)allSortLines.size() > 1)
           {
               std::sort(allSortLines.begin(), allSortLines.end(), compareNoCaseClass);
           }
        }
    }
    
    void CodBonds::getTargetBonds(std::vector<BondDict>& targetBs)
    {
        
        for (std::vector<BondDict>::iterator iCBond = allBonds.begin();
                iCBond != allBonds.end(); iCBond++)
        {
            //std::cout << iCBond->seriNum << std::endl;
            //std::cout << iCBond->atomsCodClasses[0] <<"\t" 
            //        <<  iCBond->atomsCodClasses[1] << std::endl;
            
            for (int i=0; i < (int)targetBs.size(); i++)
            {
                std::string cId0 = TrimSpaces(iCBond->atomsCodClasses[0]);
                std::string cId1 = TrimSpaces(iCBond->atomsCodClasses[1]);
                
                std::string cId2 = TrimSpaces(targetBs[i].atomsCodClasses[0]);
                std::string cId3 = TrimSpaces(targetBs[i].atomsCodClasses[1]);
                
                if(   (cId0.compare(cId2)==0 && cId1.compare(cId3)==0)
                   || (cId0.compare(cId3)==0 && cId1.compare(cId2)==0) )
                {
                    targetBs[i].codBondValues.push_back(iCBond->value);
                    //std::cout << "Bond " << iCBond->seriNum << std::endl
                    //        << "class 1" << iCBond->atomsCodClasses[0] << std::endl
                    //        << "class 2" << iCBond->atomsCodClasses[1] << std::endl
                    //        << "CCP4 dict value " << targetBs[i].length << std::endl
                    //        << "COD value " << iCBond->length << std::endl
                }
            }
        }
    }
    
    // Another class
    
    CodAngles::CodAngles()
    {
    }
    
    CodAngles::CodAngles(const CodAngles& tCAngles)
    {
        for (std::vector<AngleDict>::const_iterator iAng = tCAngles.allAngles.begin();
                iAng != tCAngles.allAngles.end(); iAng++)
        {
            allAngles.push_back(*iAng);
        }
    }
    
    CodAngles::CodAngles(FileName tFname, 
                         std::ios_base::openmode tOpenMode=std::ios::in)
    {
        if (tOpenMode == std::ios::in)
        {
            
            
            inFile.open(tFname, tOpenMode);
            
            if (inFile.is_open())
            {
                
                setupSystem();
                
                inFile.close();
                
            }
        }
    }
    
    CodAngles::~CodAngles()
    {
        if(inFile.is_open())
        {
           inFile.close(); 
        }
    }
    
    void CodAngles::setupSystem()
    {

    }
    
    void CodAngles::getTargetAngles(std::vector<AngleDict> & targetAs)
    {
       
    }
    
    
    // Another class CodTorsions
    
    CodTorsions::CodTorsions()
    {
    }
    
    CodTorsions::CodTorsions(const CodTorsions & tCT)
    {
        for (std::vector<TorsionDict>::const_iterator iTo= tCT.allTorsions.begin();
                iTo != tCT.allTorsions.end(); iTo++)
        {
            allTorsions.push_back(*iTo);
        }
    }
    
    CodTorsions::CodTorsions(FileName tFname, 
            std::ios_base::openmode tOpenMode=std::ios::in)
    {
        if (tOpenMode == std::ios::in)
        {    
            inFile.open(tFname, tOpenMode);
            
            if (inFile.is_open())
            {
                
                setupSystem();
                
                inFile.close();
                
            }
        }
    }
    
    CodTorsions::~CodTorsions()
    {
        if(inFile.is_open())
        {
            inFile.close();
        }
    }
    
    void CodTorsions::setupSystem()
    {
        std::string tRecord="";
        
        while(!inFile.eof())
        {
            std::getline(inFile, tRecord);
            tRecord = TrimSpaces(tRecord);
            std::vector<std::string> tBuf;
            StrTokenize(tRecord, tBuf);
            
            if ((int)tBuf.size() ==12)
            {
                TorsionDict aTorsion;
                
                aTorsion.seriNum  = StrToInt(TrimSpaces(tBuf[0]).substr(1,(int)tBuf[0].size()-2));
                aTorsion.value    = StrToReal(TrimSpaces(tBuf[1]));
                
                aTorsion.atomCodClasses.push_back(TrimSpaces(tBuf[2]).substr(1,(int)tBuf[2].size()-2));
                aTorsion.atomCodClasses.push_back(TrimSpaces(tBuf[3]).substr(1,(int)tBuf[3].size()-2));
                aTorsion.atomCodClasses.push_back(TrimSpaces(tBuf[4]).substr(1,(int)tBuf[4].size()-2));
                aTorsion.atomCodClasses.push_back(TrimSpaces(tBuf[5]).substr(1,(int)tBuf[5].size()-2));
                
                allTorsions.push_back(aTorsion);
                
            }
        }
    }
    
    void CodTorsions::getTargetTorsions(std::vector<TorsionDict>& targetTs)
    {
        for (std::vector<TorsionDict>::iterator iTO = allTorsions.begin();
                iTO != allTorsions.end(); iTO++)
        {
            std::cout << "COD torsion angle " << iTO->seriNum << std::endl;
            
            for (int i=0; i < (int)targetTs.size(); i++)
            {
                if ( iTO->atomCodClasses[0].compare(targetTs[i].atomCodClasses[0])==0 
                        && iTO->atomCodClasses[1].compare(targetTs[i].atomCodClasses[1])==0)
                {
                    if((iTO->atomCodClasses[2].compare(targetTs[i].atomCodClasses[2])==0 
                       && iTO->atomCodClasses[3].compare(targetTs[i].atomCodClasses[3])==0)
                       || (iTO->atomCodClasses[2].compare(targetTs[i].atomCodClasses[3])==0 
                       && iTO->atomCodClasses[3].compare(targetTs[i].atomCodClasses[2])==0) )
                    {
                        targetTs[i].codTorsionValues.push_back(iTO->value);
                        break;
                    }
                }
            }
        }
    }
}

