/* 
 * File:   plane.h
 * Author: flong
 *
 * Created on August 10, 2011, 3:40 
 */

#ifndef PLANE_H
#define	PLANE_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef RING_H
#include "ring.h"
#endif

namespace LIBMOL
{
    class Atom;
    class AtomDict;
    
    class RingDict;
    
    /* class Plane represents a plane consisting of a certain of atoms
    */
    
    class Plane   
    {
    public :
        
        // default constructor
        Plane();
        // constructor using a group of atoms
        // Plane(std::vector<Atom> tAtoms);
        // copy constructor
        Plane(const Plane & tP);
        
        // destructor
        ~Plane();

        Name getName() const;
        void setName(Name tNa);
        ID   getID() const;
        void setID(ID tID);
        SeriNumber getSeriNum() const;
        void setSeriNum(SeriNumber tSer);        
        
        Plane * create(Atom * tAtmGroup);
        void    destroy();
        
        REAL getValue(bool tN) const;
        void setValue(REAL tV, bool tS);    
                                   // the private member variables
        // void SetValue(std::vector<Atom>  tAtm);        
        REAL getSigValue() const;
        void setSigValue(REAL tV);
        REAL getValuePre() const;
        void setValuePre(REAL tV);
        
        std::vector<REAL>        atomForceConsts;
        std::vector<REAL>        planeCeofs;
        
        std::vector<Atom>        atoms;
        
    private :
        
        
  
      
        Name         itsName;
        ID           itsID;
        SeriNumber   itsSeriNum;
        
              
        
        REAL         itsValue;
        REAL         itsValueSt;
        REAL         itsSigValue;
        // REAL         itsValuePre;
        
        
    };
    
    class PlaneDict 
    {
    public :
        // Default constructor
        PlaneDict();
        
        // Copy constructor 
        PlaneDict(const PlaneDict & tPl);
        
        // Destructor 
        ~PlaneDict();
        
        ID                  archID;
        int                 archPos;
        REAL                fConst;
        std::map<ID, int>   atoms;
        
       
    };
    
    extern bool checkATorsAtomsInPl(std::vector<int> & tAtms);
    
    extern void setAllRingPlanes(std::vector<RingDict>      & tAllRings,
                                 std::vector<AtomDict >     & tAtoms, 
                                 std::vector<PlaneDict>     & tPlanes);
    
    extern void setAllRingPlanes2(std::vector<RingDict>      & tAllRings,
                                  std::vector<std::vector<int> >  & tMRingIdxs,
                                  std::vector<AtomDict >     & tAtoms, 
                                  std::vector<PlaneDict>     & tPlanes);
    
    extern void setAllOtherPlanes(std::vector<RingDict>     & tAllRings,
                                  std::vector<AtomDict >    & tAtoms, 
                                  std::vector<PlaneDict>    & tPlanes);
}




#endif	/* PLANE_H */

