/* 
 * File:   atomAssembly.h
 * Author: flong
 *
 * Created on August 9, 2011, 11:42 AM
 */

#ifndef ATOMASSEMBLY_H
#define	ATOMASSEMBLY_H

#ifndef KERNEL_H
#include "kernel.h"
#endif 

#ifndef ATOM_H
#include "atom.h"
#endif

namespace LIBMOL //temp 
{
    class Atom;
    
    /* AtomAssembly is a base class for those classes that describe the systems 
     * contain atoms, such as Bond, Angle, Model, Residue, Molecule etc.
     */
    
    class AtomAssembly : public Entity
    {
    public :
        
        // Default constructor 
        AtomAssembly();
        // Copy constructor
        AtomAssembly(const AtomAssembly & tAtomAssembly);
        // Destructor 
        virtual ~AtomAssembly();
         
        Name getName() const;
        void setName(Name tNa);
        ID getID() const;
        void setID(ID tID);
        SeriNumber getSeriNum() const;
        void setSeriNum(SeriNumber tSer);
        
        Size  getSize();      // return number of atoms in the assembly
        
        Atom & getOneAtom(int tN);
        
        // add One atom at the tPos position of the atom list
        // If the position tPos larger than the size of atoms
        // just append the tAtom to the end of the atom list
        void  addOneAtom(Atom & tAtom, unsigned int tPos);
        
        //delete the atom with SeriNumber tN 
        void  deletOneAtom(SeriNumber tN);
        
        //void  deletOneAtom(Atom & tAtm);
        
        // Replace the atom at tPos with tAtom
        void  swapOneAtom(Atom & tAtom, unsigned int tPos);
        
        bool operator == (const AtomAssembly & tAtmS) const;
        //AtomAssembly & operator =(const AtomAssembly & tAtmS);
        
        std:: vector<Atom>         atoms;
        
    private:
        
        Name                       itsName;
        ID                         itsID;
        SeriNumber                 itsSeriNum;      
        
    };
}




#endif	/* ATOMASSEMBLY_H */

