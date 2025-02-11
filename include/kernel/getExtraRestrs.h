/* 
 * File:   getExtraRestrs.h
 * Author: flong
 *
 * Created on September 26, 2011, 10:57 AM
 */

#ifndef GETEXTRARESTRS_H
#define	GETEXTRARESTRS_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef ATOMASSEMBLY_H
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

#ifndef SSBOND_H
#include "ssbond.h"
#endif

#ifndef LIBG_LINK_H
#include "libglink.h"
#endif

#ifndef CRYSTINFO_H
#include "crystInfo.h"
#endif

#ifndef SECONDARYSTRUCTURES_H
#include "secondaryStructures.h"
#endif

#ifndef PDBFILE_H
#include "PDBFile.h"
#endif

#ifndef RESTRAINTLISTS_H
#include "restraintLists.h"
#endif

#ifndef EXTRARESTRDICTFILE_H
#include "ExtraRestrDictFile.h"
#endif

#ifndef NEIGHBLIST_H
#include "neighbList.h"
#endif

namespace LIBMOL
{
    enum SearchMode
    {
        AUTO,
        USERDEF
    };
    
    enum ActType
    {
        // currently, only two of them
        ADDRESTRAINTS,
        PAIRS,
        BASEPAIRS     
    };
 
    enum AlignMode
    {
        PARA, 
        ANTIPARA
    };
    
    class MissArgException
    {
    };
    
    struct ChainSect
    {
        int  from;
        int  to;
        ID   chainID;
    };
    
    class restrListsForOnePair
    {
    public :
        restrListsForOnePair();
        restrListsForOnePair(const restrListsForOnePair & tRLFOP);
        
        // sequence for sortKeys [0] chainID, [1] seqNum, other sort keys 
        // could be added into the vector when needed 
        std::vector<int> sortKeys;   
        
        std::string pairKey;
        std::vector<ResidueID>       resPair;
        std::vector<RestraitListSet> restrLists;
        
    };
    
    bool SortByKeys(const restrListsForOnePair  & obj1,
                    const restrListsForOnePair  & obj2);
    
    class userDefActSet
    {
    public :
        // default constructor
        userDefActSet();
        // copy constructor 
        userDefActSet(const userDefActSet & tU);
        // destructor
        ~userDefActSet();
        
        SearchMode  searchMode;
        ActType     actType;
        std::vector<ChainSect> residueSet;
        AlignMode   alignMode;
        std::string label;
    };
    
    class getUserDefAction 
    {
    public :
        
        // default constructor 
        getUserDefAction();
        
        // constructor using a file name
        getUserDefAction(FileName  uF);
        
        // destructor 
        ~getUserDefAction();
        
        void getAllUserDef(FileName uF);
        
        std::vector<userDefActSet>                   userDefActs;       
        std::vector<std::pair<Residue, Residue> >    residuePairs;
    private: 
        
    };   
    
    class GetExtraRestrs 
    {
    public :
        
        // default constructor 
        GetExtraRestrs();
        // copy constructor
        GetExtraRestrs(const GetExtraRestrs & tG);
        // constructor with three file name inputs
        // GetExtraRestrs(Name inPDBName, Name inDictName, 
        //               Name outFileName, Name inMonoList, Name uDefF);
 
        // constructor with three file name inputs of c-string (char *)
        GetExtraRestrs(FileName inPDBName, FileName inDictName, 
                       FileName outFileName, FileName inMonoList,
                       FileName uDefF);
        
        // destructor
        ~GetExtraRestrs();
        
        void       execute();
        
        int        getDim();
        void       setDim(int tN);
        
        SearchMode getSearchMode();    
        void       setSearchMode(SearchMode tS);
        ActType    getActType();
        void       setActType(ActType tA);
        
        void setResidueGroupType();
        bool checkReducedGroups(std::vector<Residue>::iterator iR);
        void checkFullGroups(std::vector<Residue>::iterator iR);
        
        void autoFindExtraRestraints();
        void autoFindBasePairList(NeighbList & tNBList);
        void setupDNARNAList();
        bool checkDNARNABase(Name tNa);
        bool checkDNARNABasePair(Name tR1, Name tR2);
        
        void findExtraRestraints();
        void searchPairRestraints(int i1, ID tID1, int i2, ID tID2);
        Residue & findOneResidue(int tN, ID iID);
        void findOneResidue(int tN, ID iID,
                            std::vector<Residue> & tPair);
        void buidResidueSet(ChainSect & tCh, 
                            std::vector<Residue> & tResSet);
        
        void findAllRestraintsForOnePair(std::vector<Residue> & tPair);
        void findAllAddRestraintsForOneResidue(
                 std::vector<Residue>::iterator  iR,
                 std::string  tS);
        
        void searchTorsions();
        void searchAddTorsions(std::vector<Residue>::iterator  tR,
                               std::vector<Torsion>::iterator iR,
                               Name tLabel);
        void searchOnePairTorsion(std::vector<Residue>  & tR,
                               std::vector<Torsion>::iterator iR,
                               Name tLabel);
        
        
        void searchBonds();
        void searchAddBonds(std::vector<Residue>::iterator tR,
                            std::vector<Bond>::iterator iB,
                            Name tLabel);
        void searchOnePairBond(std::vector<Residue>  & tR,
                               std::vector<Bond>::iterator iB,
                               Name tLabel);
        
        void searchAngles();
        void searchAddAngles(int tN);     
        void searchOnePairAngle(std::vector<Residue>  & tR,
                                std::vector<Angle>::iterator iA,
                                Name tLabel);
        
        
        void searchChirals();
        void searchAddChirals(int tN); 
        void searchOnePairChiral(std::vector<Residue>  & tR,
                                std::vector<Chiral>::iterator iA,
                                Name tLabel);        
        
        void searchPlanes();
        void searchAddPlanes(int tN);
        void searchPairPlanes(int tN);
        void searchOnePairPlane(std::vector<Residue>  & tR,
                                std::vector<Plane>::iterator iP,
                                Name tLabel);
        
        void makeResiduePairs();
        
        
        std::string   getErrInfo() const;
        void          setErrInfo(std::string tS);
        
        int           getErrLevel() const;
        void          setErrLevel(int tN);
        
        // check and validate the restraints found
        void checkAllRestrFound();
        void checkRestraintsFromOnePair(
               std::map<std::string, RestraitListSet>::iterator iM);
        void checkRestrListSize(
               std::map<std::string, RestraitListSet>::iterator iM);
        void copyRestraintPair(
               std::map<std::string, RestraitListSet>::iterator iM);
        void reOrderAllRestrFound();
        
        
        // output the restraints found 
        void outputExtraRestraints(int tMode);
        void outputAddRestraints();
        void outputBondRestraints(std::ofstream & tO,
                                  std::vector<Bond> & tBList,
                                  int tMode);
        void outputAngleRestraints(std::ofstream & t,
                                   std::vector<Angle> & tAList,
                                   int tMode);
        void outputTorsionRestraints(std::ofstream & tO,
                                     std::vector<Torsion> & tTList,
                                     int tMode);
        void outputChiralRestraints(std::ofstream & tO,
                                    std::vector<Chiral> & tRList,
                                    int tMode);
        
        void outInstructFile();
        
        // void outputPlanneRestraints(std::ofstream & tO);
        // input data 
        PDBFile            *    PDBObj;
        ExtraRestrDictFile *    restrDictObj; 
        MonomerLib         *    monomerList;
        getUserDefAction   *    allUserDefAct;
        std::vector<std::pair<Residue, Residue> >   residuePairs;
        std::vector<NeighbList>  allNBLists;
        // neighbor lists
        
        std::map<std::string, std::string>    reducedGroups;
        std::map<std::string, RestraitListSet> allfoundRestrLists;
        std::map<std::string, RestraitListSet> allProbRestrLists;
        std::vector<restrListsForOnePair> sortedAllRestLists;
        FileName              outRestrFileName;

        // some data 
        std::vector<std::string>                  DnaRnaBaseList;
        std::multimap<std::string, std::string>   DnaRnaBasePairs;
        // some exception classes 
               // exception classes 
        class findRestrException
        {
        };
        class outRestrException
        {
        };
        class noMonoGroupException : public std::exception
        {
            virtual const char* what() const throw()
            {
                return "No group type is found for the residue";
            }
        };
        
    private :
        
        int                     itsDim;
        SearchMode              itsSearchMode;
        ActType                 itsActType;
        std::string             itsErrInfo;
        int                     itsErrLevel;   
        
    };
    

}

#endif	/* GETEXTRARESTRS_H */

