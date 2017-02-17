/* 
 * File:   neighbList.h
 * Author: flong
 *
 * This code builds cell-linked lists and neighbor lists for every atom
 * in the system.
 * Created on August 11, 2011, 5:54 PM
 */

#ifndef NEIGHBLIST_H
#define	NEIGHBLIST_H

#ifndef KERNEL_H
#include "kernel.h"
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

#ifndef CRYSTINFO_H
#include "crystInfo.h"
#endif

namespace LIBMOL
{
    class Atom;
    class AtomDict;
    
    // class NeighbList describes a list of atoms neighboring to a particle
    // atom
    
     struct ResidueID
    {
        Name          resName;
        SeriNumber    seqNum;
        ID            chainID;
    };   
 
    class NBCell
    {
    public :
        NBCell();
        NBCell(const NBCell & tC);
        ~NBCell()
        {
        }
        
        std::vector<int>                 index;
        std::vector<Atom>                atomsInCell;
        std::vector<ResidueID>           residueList;
        std::vector<std::vector<int> >   nbcellList;
    };
    
  
    class NeighbList 
    {
    public :
        
        // default constructor
        NeighbList();
        // copy constructor
        NeighbList(const NeighbList & tNBL);
        // constructor with detailed information
        NeighbList(std::vector<Atom> & tAtms,
        int tDim, REAL tNBCutoff, REAL tNBShell);
        
        //destructor 
        ~NeighbList();
        
        // exception classes
        class cellLengthException
        {
        };
        
        
        
        int         getDim() const;
        void        setDim(int tD);
        
        REAL        getCutOff() const;
        void        setCutOff(REAL tC);
        
        REAL        getNBShell() const;
        void        setNBShell(REAL tC);
                
        std::string getErrInfo() const;
        void        setErrInfo(std::string tE);
        
        int         getErrLevel() const;
        void        setErrLevel(int tN);
                
        void        buildCellNBList(std::vector<NBCell>::iterator tNB);
        void        putAtomsAndResiduesInCell();
        void        buildCellSystem();
        
        void        buildLinkedList();
        SeriNumber  getCellID(SeriNumber tIdx); 
        int         getOneNBCellPos(std::vector<int> & tIdx);
        bool        getOneNBCell(std::vector<int> & tIdx);
        bool        getOneNENBCell(std::vector<int> & tIdx);
        void        buildAtomNeighbList();
        void        buildAtomNeighListfromAPairCell(NBCell & tCell1,
                                                    NBCell & tCell2);
        
        void        updateResidueNBList(std::vector<Atom>::iterator tA);
        void        buildResidueNBList();
        void        buildResidueNeighbListfromAPairCell(NBCell & tCell1,
                                                        NBCell & tCell2);
        void        building(int tMode);   
        
        std::vector<Atom>                            allAtoms;
        std::map<int, std::list<Atom> >              atomNBlist;
        std::map<std::string, std::list<ResidueID> > residueNBList;   
        std::vector<REAL>                            coordsMax;
        std::vector<REAL>                            coordsMin;
        std::vector<NBCell>                          allCells;
        std::vector<NBCell>                          allNECells;
        int                                          totalNumCells;
        
    private :
        
        // parameters
        int                 itsDim;
        REAL                itsCutOff;
        REAL                itsNBShell;
        
        std::vector<REAL>   itsCellLength;
        std::vector<int>    itsNumCell;
        std::string         itsErrInfo;
        int                 itsErrLevel;
             
    };
    
    class NBCellDict
    {
    public :
        NBCellDict();
        NBCellDict(const NBCellDict & tC);
        ~NBCellDict()
        {
        }
        
        std::vector<int>                 index;
        std::vector<int>                 atomsInCell;
        std::vector<ResidueID>           residueList;
        std::vector<std::vector<int> >   nbcellList;
    };
    
    class NeighbListDict 
    {
    public :
        
        // default constructor
        NeighbListDict();
        
        //destructor 
        ~NeighbListDict();
        
        // exception classes
        class cellLengthException
        {
        };
        
        
        void        buildCellNBList(std::vector<NBCellDict>::iterator tNB);
        void        putAtomsAndResiduesInCell(std::vector<AtomDict> & aAtomList);
        void        buildCellSystem(std::vector<AtomDict> & aAtomList,
                                    int tDim, REAL tNBCutoff, REAL tNBShell);
        void        buildAtomNeighbList(std::vector<AtomDict> & tAtomList, REAL tL);
        void        buildAtomNeighListfromAPairCell(std::vector<AtomDict> & tAtomList,
                                                    REAL tL,
                                                    NBCellDict & tCell1,
                                                    NBCellDict & tCell2);
        
        void        updateResidueNBList(std::vector<AtomDict>::iterator tA);
        void        buildResidueNBList();
        void        buildResidueNeighbListfromAPairCell(NBCellDict & tCell1,
                                                        NBCellDict & tCell2);
        void        building(std::vector<AtomDict> & aAtomList,
                             int tDim, REAL tNBCutoff, 
                             REAL tNBShell, int tMode);  
        void        building(std::vector<AtomDict>   &   tAtomList, 
                             std::vector<AtomDict>   &   tAllAtoms,
                             int tDim,              REAL tNBCutoff, 
                             REAL tNBShell,              int tMode);
        
        
       
        std::map<std::string, std::list<ResidueID> > residueNBList;   
        std::vector<REAL>                            coordsMax;
        std::vector<REAL>                            coordsMin;
        std::vector<NBCellDict>                      allCells;
        std::vector<NBCellDict>                      allNECells;
        int                                          totalNumCells;
        
        
        // parameters
        
        
        std::vector<REAL>   itsCellLength;
        std::vector<int>    itsNumCell;
        std::string         itsErrInfo;
        int                 itsErrLevel;   
        
    };
}

#endif	/* NEIGHBLIST_H */

