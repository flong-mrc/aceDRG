/*
 * File:   chain.h
 * Author: flong
 *
 * Created on August 11, 2011, 1:34 PM
 */

#ifndef CHAIN_H
#define	CHAIN_H

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

namespace LIBMOL
{
    class Atom;
    class Residue;
    
    /* Class chain describes a polypeptide chain in a protein model 
     * 
     */
    class Chain 
    {
    public :
        
        //default constructor
        Chain();
        // copy constructor
        Chain(const Chain & tC);
     
        // constructor with mini info
        Chain(ID tID, SeriNumber tN, Name tName, std::vector<Residue> & tResG);
        
        // destructor
        virtual ~Chain();

        Chain & operator =(const Chain & tC);
        

        Name getName() const;
        void setName(Name tNa);
        ID getID() const;
        void setID(ID tID);
        SeriNumber getSeriNum() const;
        void setSeriNum(SeriNumber tSer);
        SeriNumber getModSeriNum() const;
        void  setModSeriNum(SeriNumber tN);
        
        Atom *      getCurrentAtom() const;
        void        setCurrentAtom(const Atom &  tAtom);
        Size        getNumOfRes() const;
        void        setNumOfRes(int iN); 
        
        Residue  &  getOneResidue(SeriNumber tSeqNum);
        
        // add a residue at the end of the residue list
        void        addOneResidue(Residue & tR);
        
        // add a residue at the position tPos
        void        addOneResidue(Residue & tR, unsigned int tPos); 
        
        // void        addOneResidue(ID tID, SeriNumber tN, Name tName);
        
        // delete one residue from the chain at the position tPos
        void        deleteOneResidue(unsigned int tPos);
        // delete the residue with tID, tN, and tName from the chain 
        void        deleteOneResidue(ID tID, SeriNumber tN, Name tName);
        
        // The tN-th residue in the chain is replaced with residue tR
        void        swapOneResidue(unsigned int tN, Residue & tR);
        
        
        bool operator==(const Chain& chain) const;
    
        std::map<std::string, bool>     state;
        /* key-type for state are:
         * isAminoAcid
         * isNucleotide
         * isSolvent
        */
        
        std::vector<Residue>       residues;
        SeriNumber                 curAtomSeriNum;
    private:

        /* defined in class AtomAssembly*/
        Name              itsName;
        ID                itsID;
        SeriNumber        itsSeriNum;
        
        int               itsNumOfResidues;
        
        SeriNumber        itsModSeriNum;
        
        Atom          *   itsCurrentAtom;
        
        
    };

}

#endif	/* CHAIN_H */

