/* 
 * File:   Residue.h
 * Author: flong
 *
 * Created on August 10, 2011, 5:02 PM
 */

#ifndef RESIDUE_H
#define	RESIDUE_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef LIBG_MODEL_H
#include "libgmodel.h"
#endif

#ifndef LIBG_ATOMASSEMBLY_H
#include "atomAssembly.h"
#endif

namespace LIBMOL
{
    class Atom;
    class AtomDict;
    class Chain;
    class Model;
    class Molecule;
    class Fragment;
    
    class TorsionList;
    class AngleList;
    class BondList;
    
    /* Class Residue represent a specific monomer within the polymeric chains.
     * Here we mainly refer to amino acid residues and residues in RNA/DNA
     */
    
    class Residue   
    {
    public :
        
        // default constructor 
        Residue();
        // copy constructor
        Residue(const Residue & tRes);
        // constructor with mini info
        Residue(ID tID, SeriNumber tN, Name tName);
        
        // destructor
        ~Residue();
        
       // Residue & operator=(const Residue & tR);
        // create a new residue
        Residue & create(ID tID, SeriNumber tN, Name tName, Atom * tAtmGrp);
        
        // add, delete and sway atom etc have taken 
        // care of by class atomAssembly such as 
        //   Atom *          getAtomGroup();
          
        Name getName() const;
        void setName(Name tNa);
        ID getID() const;
        void setID(ID tID);
        SeriNumber getSeriNum() const;
        void setSeriNum(SeriNumber tSer);
        
        SeriNumber  getSeqNum() const;
        void setSeqNum(SeriNumber tN);
        
        ID getInsCode() const;
        void setInsCode(ID tIn);
        
        SeriNumber getChainSeriNum() const;
        void   setChainSeriNum(SeriNumber tN);
 
        ID getChainID() const;
        void setChainID(ID tID);

        SeriNumber  getModelSeriNum() const;
        void        setModelSeriNum(SeriNumber tN);
        
    //    Atom  & getCurrentAtom() const;
    //    void    setCurrentAtom(Atom & tAtom);
 
        ID   getGroupID() const;
        void setGroupID(ID tID);

        void addOneAtom(Atom &tAtom);
        
        
              
        // public member variables
        std::map<std::string, bool>     state;
        /* key-type for state are:
         * isAminoAcid
         * isCTerminal
         * isNTerminal
         * isNucleotide
         * isDNARNA
         * isSugar
         * isSolvent
         */
        
        std::vector<Atom>           atoms;
       // std::vector<Atom>  itsCurrentAtom;
        
    private : 
        
        /* defined in class AtomAssembly*/
        Name              itsName;
        ID                itsID;
        SeriNumber        itsSeriNum;
        
        
        SeriNumber        itsSeqNum;
        ID                itsInsCode;
        SeriNumber        itsChainSeriNum;
        ID                itsChainID;
        SeriNumber        itsModelSeriNum; 
        
        ID                itsGroupID;
        
        
    };
    
    class ModRes : public Residue 
    {
    public:
        
        // default constructor 
        ModRes();
        
        //copy constructor
        ModRes(const ModRes  & tM);
        
        // destructor
        ~ModRes();
        
        
        inline Name getStdName() const
        {
            return itsStdName;
        }
        inline void setStdName(Name tName)
        {
            itsStdName = tName;
        }
        
        inline std::string getComment() const
        {
            return itsComment;
        }
        inline void setComment(std::string tC)
        {
            itsComment = tC;
        }
        
        ID getInsCode() const;
        void setInsCode(ID tS);   
        // std::vector<Atom>    itsCurrentAtom;
   private:
        Name              itsName;
        ID                itsID;
        SeriNumber        itsSeriNum;
        
        
        SeriNumber        itsSeqNum;
        ID                itsInsCode;
        SeriNumber        itsChainSeriNum;
        ID                itsChainID;
        SeriNumber        itsModelSeriNum; 
        
        Name              itsStdName;
        std::string       itsComment;
        
        
        
    };
    
    class ResidueDict 
    {
    public:
        // default constructor 
        ResidueDict();
        // copy constructor
        ResidueDict(const ResidueDict & tRes);
        // constructor with mini info
        ResidueDict(SeriNumber tN, Name tName, 
                    SeriNumber tSeqNum, ID tChainID);
        
        // destructor
        ~ResidueDict();
        
        SeriNumber                  seriNum;
        Name                        name;
        SeriNumber                  seqNum;
        ID                          chainID;
        
        std::vector<AtomDict>       atoms;
        
    };
    
}
    



#endif	/* RESIDUE_H */

