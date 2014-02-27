/* 
 * File:   Bond.h
 * Author: flong
 *
 * Created on August 9, 2011, 11:31 AM
 */

#ifndef BOND_H
#define	BOND_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef LIB_ATOMASSEMBLY_H
#include "atomAssembly.h"
#endif

#ifndef RESIDUE_H
#include "residue.h"
#endif

namespace LIBMOL
{
    class Atom;
    class AtomDict;
    class Residue;
    /* Bond class describes bonding between two or more atoms. It derived from
     * the base class AtomAssembly
     * Bond class has the following components, besides constructors and
     * destructor,  
     * 
     * @ member variables   
     * 
     */ 
    
    class Bond 
    {
    public :
        

        typedef short       Order;
        typedef short       Type;
        
        // default constructor 
        Bond();
        
        // copy constructor 
        Bond(const Bond & tBond);
        
        // create bonding constructor
        Bond(Atom & tAtm1, Atom & tAtm2);
        
        // for disulfide bond
        Bond(Residue & tRes1, Residue & tRes2);
        
        // destructor
        virtual ~Bond();
        
        Bond & operator = (const Bond & tB);
   
        Name getName() const;
        void setName(Name tNa);
        ID getID() const;
        void setID(ID tID);
        SeriNumber getSeriNum() const;
        void setSeriNum(SeriNumber tSer);
        
        //Bond * createBond(Atom& tAtm1, Atom& tAtm2);
        Bond * createBond(Atom & tAtm);
        void   swapAtom(Atom &tAtm);
        
        void   destroyBond();   // set everything to null
        
        ID    getOrder(bool tSt) const;
        void  setOrder(ID tO, bool  tSt);

        
        // this function set the bond-order according to private 
        // variable 'itsStates'
        
        //void setOrder();
        
        inline Type getType () const
        {
            return itsType;
        }
        inline void setType(Type tT)
        {
            itsType = tT;
        }
        // this function set the bond-order according to private 
        // variable 'itsStates'
        //void setType();
        REAL getLength(bool tSt) const;
        void setLength();
        void setLength(Atom & tAtm1, Atom & tAtm2);
        void setLength(REAL tL, bool tSt);
        
        inline REAL getSigLength() const
        {
            return itsSigLength;
        }
        inline void setSigLength(REAL tS)
        {
            itsSigLength = tS;
        }
        
        inline REAL getForceConst() const 
        {
            return itsForceConst;
        }
        inline void setForceConst(REAL tF)
        {
            itsForceConst = tF;
        }
        
        
        REAL getAtomPairLength(Atom & tA1, Atom tA2);
        
        std::vector<Atom>             atoms;  
        bool                          isItTouched;
        std::map<std::string, bool>   states;    
        
        
        
    private :
        
        // Assume two atoms construct a bond
           
        
        /* defined in class AtomAssembly */
        Name           itsName;
        ID             itsID;
        SeriNumber     itsSeriNum;
        
        
        ID             itsOrder;
        ID             itsOrderSt;
                
        Type           itsType;
        
        // the following map decides what order and type of the bonding is       
        
        REAL           itsLength;
        REAL           itsLengthSt;
        REAL           itsSigLength;     
        REAL           itsForceConst;     
        
    };
    
    class BondDict 
    {
    public:
        
        // Default constructor 
        BondDict();
        // copy constructor 
        BondDict(const BondDict & tBond);
        // destructor 
        ~BondDict();
                
        std::string bondOrderNumToStr();
        
        ID                      resName;
        SeriNumber              seriNum;
        
        ID                      order;
        REAL                    value;
        REAL                    sigValue;
        REAL                    valueST;
        REAL                    sigValueST;
        REAL                    valueP;
        REAL                    sigValueP;
        
        bool                    hasMetal;
        bool                    hasCodValue;
        int                     numCodValues;
        int                     numCodValuesP;
        
        std::map<ID, int>       fullAtoms; // for atoms from SMILE, ID is chemType
        std::vector<ID>         atoms;
        std::vector<int>        atomsIdx;
        
        std::vector<ID>         atomsCodClasses;
        std::vector<int>        atomsHashingCodes;
        
        std::vector<ID>         atomsNBRep;
        std::vector<ID>         atomsNB2Rep;
        std::vector<REAL>       codBondValues;
       
        
    };
    
    class MetBond 
    {
    public:
                   // Default constructor 
        MetBond();
        // copy constructor 
        MetBond(const MetBond & tBond);
        // destructor 
        ~MetBond();           
        
        ID                      resName;
        SeriNumber              seriNum;
        
        short                   order;
        REAL                    value;
        REAL                    sigValue;
        REAL                    valueST;
        REAL                    sigValueST;
        
        
        bool                    hasCodValue;
        int                     numCodValues;
        
        
        int                     mAtomIdx;
        int                     lAtomIdx;
        
    };
    
    extern int getBond(std::vector<BondDict> & tAllBonds, int tAt1, int tAt2);
}



#endif	/* BOND_H */

