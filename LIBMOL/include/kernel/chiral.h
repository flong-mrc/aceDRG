/* 
 * File:   chiral.h
 * Author: flong
 *
 * Created on August 10, 2011, 2:14 PM
 */

#ifndef CHIRAL_H
#define	CHIRAL_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef ATOMASSEMBLY_H
#include "atomAssembly.h"
#endif

namespace LIBMOL
{
    class Atom;
    
    /* Class Chiral represents chirality in molecular geometry 
     *  
    */
    
    class Chiral 
    {
        
    public :
        
        // default constructor
        Chiral();
        // copy constructor 
        // Chiral(Chiral & tC);
        
        Chiral(const Chiral & tC);
        
        
        // destructor 
        ~Chiral();
        
        //Chiral & operator=(const Chiral &tC);
 
        Name getName() const;
        void setName(Name tNa);
        ID getID() const;
        void setID(ID tID);
        SeriNumber getSeriNum() const;
        void setSeriNum(SeriNumber tSer);
        Chiral create (Atom*  atmList);
        void   destroy();
 
        REAL getValue(bool tL) const;
        void setValue(REAL tV, bool tL);
        void setValue();             
        REAL setValue(Atom & tA1, Atom & tA2, Atom & tA3, Atom & tA4);
        
        // void SetValue(std::vector<Atom> tAtms);        
        REAL getSigValue() const;
        void setSigValue(REAL tS);
        REAL getValuePre() const;
        void setValuePre(REAL tVP);
        
        REAL getForceConst() const;
        void setForceConst(REAL tF);
        
        bool                    isItTouched;
        std::vector<Atom>       atoms;
   
        
    private:
        
        Name           itsName;
        ID             itsID;
        SeriNumber     itsSeriNum; 
        
        
        REAL                    itsValue;
        REAL                    itsValueSt;
        REAL                    itsSigValue;
        REAL                    itsValuePre;
     
        REAL                    itsForceConst;
        
    };
    
    // Another chiral class for the dictionary
    class ChiralDict
    {
    public:
        
        // Default constructor
        ChiralDict();
        
        // Copy constructor
        ChiralDict(const ChiralDict & tCh);
        
        // Destructor
        ~ChiralDict();
        
        ID                  id;       
        ID                  archID;
        int                 archPos;
        REAL                value;
        ID                  sign;
        REAL                valueST;
        ID                  signST;
        REAL                fConst;
        //std::map<ID, int>   atoms;
        std::vector<int>    atoms;
        std::map<int, std::vector<int> > mutTable;
        
    };
    
    extern int inChirals(std::vector<ChiralDict> tChirals, 
                         AtomDict & tInAtom);
    
    extern void buildChiralCluster(std::vector<int>  & tInAts,
                                   std::vector<int>  & tOutAts,
                                   int                 tExcluded,
                                   int                 tChstat);
    
    extern void buildChiralCluster2(ChiralDict        & tChiral,
                                    std::vector<int>  & tOutAts, 
                                    int                 tExcluded);
    extern void getMutList(std::vector<int> & tIn, std::vector<int> & tOut);
}


#endif	/* CHIRAL_H */

