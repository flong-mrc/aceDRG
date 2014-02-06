/* 
 * File:   Torsion.h
 * Author: flong
 *
 * Created on August 9, 2011, 11:32 AM
 */

#ifndef TORSION_H
#define	TORSION_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

namespace LIBMOL
{
    
    class Atom;
    class AtomDict;
    
    /* class Torsion describes torsion angles formed by 4 atoms
     *
     * 
     *
     */
    class Torsion 
    {
    public:
        
        // default constructor
        Torsion();
        // copy constructor
        Torsion(const Torsion & tT);
        // constructor by 4 atoms
        Torsion(Atom& tAtm1, Atom & tAtm2, Atom& tAtm3, Atom& tAtm4);
        
        // destructor
        virtual ~Torsion();

        Name getName() const;
        void setName(Name tNa);
        ID getID() const;
        void setID(ID tID);
        SeriNumber getSeriNum() const;
        void setSeriNum(SeriNumber tSer);
        
        Torsion create(Atom& tAtm1, Atom& tAtm2, Atom& tAtm3, Atom& tAtm4);
        void    destroy();
 
        REAL getValue(bool tB) const;
        void setValue(REAL tT, bool tB); 
        void setValue();       
        // calculate the torsion angle value
        // from 4 atoms       
        REAL setValue(Atom & tA1, Atom & tA2, Atom & tA3, Atom & tA4);  
        
        REAL getSigValue() const;
        void setSigValue(REAL tS);
        
        REAL getValuePre() const;
        void setValuePre(REAL tVP);
        
        REAL getForceConst() const;
        void setForceConst(REAL tF);
        
        int  getPeriod() const;
        void setPeriod(int tP);
        
 
  
        // here are a group variables without initialization 
        // when the object is created
        REAL  disAtmSE;
        REAL  disAtmSEDev;
        REAL  velo;
        REAL  veloPre;
        REAL  accel;
        REAL  accelPre;
        
        
        // Variables linked only with a tree structure
        
        int        varFlag;  
                                     
        bool       isItTouched;
        
        std::vector<REAL>     normalUnit;
        std::vector<int>      torsUpAtomList;
        std::vector<Atom>     atoms;
        
        std::map<std::string, bool> states;
        /* main key-types are
         * isPhi
         * isPsi
         * isOmg    
         */  
        
        
        
    private:
        
       
        Name           itsName;
        ID             itsID;
        SeriNumber     itsSeriNum;
           
        
        REAL         itsValue;
        REAL         itsValueSt;
        REAL         itsSigValue;
        REAL         itsValuePre;
        
        REAL         itsForceConst;    
        
        int          itsPeriod;  
        
        
               
    };
    
    class TorsionDict
    {
    public :
        
        //Default constructor
        TorsionDict();
        // Copy constructor 
        TorsionDict(const TorsionDict & tTorsion);
        // Constructor using atoms
        TorsionDict(std::vector<AtomDict> & tAtoms);
        // Default destructor
        ~TorsionDict();
        
        int                   seriNum;
        REAL                  value;
        REAL                  sigValue;
        REAL                  valueST;
        REAL                  sigValueST;
        int                   period;
        ID                    id;
        
        std::vector<int>      atoms;
        std::vector<ID>       atomCodClasses;
        std::vector<REAL>     codTorsionValues;  
        
    };
    
    // Get the index of a torsion in a vector of torsions
    extern int getTorsion(std::vector<TorsionDict> & tTors,
                          int tAt1, int tAt2, int tAt3, int tAt4);
    
    // Get the value of the existing torsion based atom tree 
    extern REAL getTorsion(std::vector<AtomDict> & tAtoms,
                           int iCur, int iNext, std::vector<int> tDoneSet);
    // calculate the value of a torsion base on coords of its 4 atoms
    extern REAL getTorsion(AtomDict & tA1, 
                           AtomDict & tA2, 
                           AtomDict & tA3, 
                           AtomDict & tA4);
}

#endif	/* TORSION_H */

