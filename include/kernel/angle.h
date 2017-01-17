/* 
 * File:   Angle.h
 * Author: flong
 *
 * Created on August 9, 2011, 11:32 AM
 */

#ifndef ANGLE_H
#define	ANGLE_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef ATOMASSEMBLY_H
#include "atomAssembly.h"
#endif

#ifndef BOND_H
#include "bond.h"
#endif

#ifndef UTILITY_H
#include "utility.h"
#endif

namespace LIBMOL // temp 
{
    class Atom;
    class AtomDict;
    class Bond;
    
    /* Angle class describes the angle formed among three atoms *
     * 
     * 
     */
    
    class Angle
    {
    public:
        
        // default constructor
        Angle();
        // copy constructor 
        Angle(const Angle& tA);
        // constructor using 3 atoms
        Angle(Atom& tAtm1, Atom& tAtm2, Atom tAtm3);
        
        // destructor
        ~Angle();
        
        Angle & operator = (const Angle & tA);

        Name getName() const;
        void setName(Name tNa);
        ID   getID() const;
        void setID(ID tID);
        SeriNumber getSeriNum() const;
        void setSeriNum(SeriNumber tSer);
        
        Angle & create(Atom& tAtm1, Atom& tAtm2, Atom& tAtm3);
        void    swapOneAtom(Atom & tAtm);
        void    destroy();   // set everything to null
        

        REAL getValue(bool tL) const;
        void setValue();
        void setValue(REAL tV, bool tL);
        void setValue(Atom& tA1, Atom& tA2, Atom tA3);
        
        REAL getSigValue(bool tL) const;
        void setSigValue(REAL tV, bool tL);
        
        inline REAL getForceConst() const
        {
            return itsFrceConst;
        }
        inline void setForceConst(REAL tF)
        {
            itsFrceConst = tF;
        }
       
        REAL getValuePre() const;
        void setValuePre(REAL tVP);
       
        
        // set those member variables public for efficiency  
        
        bool            isItTouched;
        
        
        std::vector<Atom>    atoms;
        
    private:
        
        
        
        
        Name           itsName;
        ID             itsID;
        SeriNumber     itsSeriNum;     
        
        REAL            itsValue;
        REAL            itsValueSt;
        REAL            itsSigValue;
        REAL            itsSigValueSt;
        
        REAL            itsFrceConst;           
        
    };
    
    class AngleDict 
    {
    public:
        // Default constructor 
        AngleDict();
        // Copy constructor 
        AngleDict(const AngleDict & tAngle);
        // Constructor using atoms
        AngleDict(ID tAnchID, int tAnchPos, std::vector<AtomDict> & tAtoms);
        // Default destructor 
        ~AngleDict();
        
        // void setValue();
        
        
        
        int                     seriNum;
        REAL                    value;
        REAL                    sigValue;
        REAL                    valueST;
        REAL                    sigValueST;
        bool                    hasCodValue;
        int                     numCodValues;
        REAL                    valueP;
        REAL                    sigValueP;
        int                     numCodValuesP;
        
        ID                      anchorID;
        int                     anchorPos;
        
        bool                    isFixed;
        int                     approxLevel;
        int                     isInSameRing;
        
        std::vector<int>        atoms;
        std::vector<ID>         atomChemTypes;
        std::vector<ID>         atomsCodClasses;
        std::vector<ID>         atomsNB2Rep;
        std::vector<ID>         atomsNBRep;
        std::map<ID, ID>        atomsSPStats;
        std::map<ID, ID>        atomsNB1NB2SPStats;
        std::vector<REAL>       codAngleValues;
    };
    
    extern int getAngle(std::vector<AngleDict> & tAllAngles,
                        int cAt, int tAt1, int tAt2);
    
    extern REAL getAngleValueFromFracCoords(AtomDict  & tAtCen,
                                            AtomDict  & tAt1, 
                                            AtomDict  & tAt2,
                                            REAL a, REAL b, REAL c, 
                                            REAL alpha, REAL beta, REAL gamma);
    
}

#endif	/* ANGLE_H */

