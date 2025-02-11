/* 
 * File:   secondaryStructures.h
 * Author: flong
 *
 * Created on September 7, 2011, 6:52 PM
 */

#ifndef SECONDARYSTRUCTURES_H
#define	SECONDARYSTRUCTURES_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef FILE_H
#include "file.h"
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

namespace LIBMOL
{
    class Atom;
    class Residue;
    class Chain;
    class Model;

    enum TypeOfHelix
    {
        Unknown            =   0,
        Right_handed_alpha =   1,
        Right_handed_omega,
        Right_handed_pi, 
        Right_handed_gamma,
        Right_handed_3_10, 
        Left_handed_alpha,
        Left_handed_omega,
        Left_handed_gamma,
        ribbon_helix_2_7,
        Polyproline,
    };

    class Helix
    {
        
    public :
        
        // default constructor
        Helix();
        
        // copy constructor
        Helix(const Helix & tHelix);
        
 
        // destructor 
        ~Helix();
       
        inline SeriNumber getSeriNum()
        {
            return itsSeriNum;
        }
        inline void setSeriNum(SeriNumber tN)
        {
            itsSeriNum = tN;
        }
        
        inline ID getID()
        {
            return itsID;
        }
        inline void setID(ID tID)
        {
            itsID = tID;
        }
        
        inline int  getHelixClass()
        {
            return  itsHelixClass;
        }
        
        inline void  setHelixClass(int tH)
        {
            itsHelixClass = tH;
        }
        
        inline std::string getComment()
        {
            return itsComment;
        }
        inline void setComment(std::string tC)
        {
            itsComment = tC;
        }
        
        inline int getLength()
        {
            return itsLength;
        }
        inline void setLength(int tL)
        {
            itsLength = tL;
        }
        
        
        Residue  *   getInitResidue();
        void         setInitResidue(Residue & tRes);
        Residue  *   getEndResidue();
        void         setEndResidue(Residue & tRes);
        
        std::vector<Residue>      residues;
        SeriNumber                itsSeriNum; 
        
 // private :
        
        
        ID            itsID;
        
        Residue   *   itsInitRes;
        Residue   *   itsEndRes;
        
        int           itsHelixClass;
        std::string   itsComment;
        
        int           itsLength;
                     
        
    };
    
    
    class Strand 
    {
    public :
        
        // default constructor
        
        Strand();
        
        // copy constructor 
        Strand(const Strand & tStrand);
        
        // destructor 
        ~Strand();
        
        inline SeriNumber getSeriNum()
        {
            return itsSheetSeriNum; 
        }
        inline void setSeriNum(SeriNumber tN)
        {
            itsSheetSeriNum = tN;
        }
        
        inline ID getID()
        {
            return itsID;
        }
        inline void setID(ID tID)
        {
            itsID = tID;
        }
              
        // Get or set the initial or ending residues. 
        // 1 for initial one and -1 for the ending one    
        Residue  *   getTerResidue(int tN);
        void         setTerResidue(Residue & tRes, int tN);
    
       // Get or set the current or previous registration atoms in the strand
       // 1 for the current one and -1 for the previous one
       Atom      *   getRegistAtom(int tN);
       void         setRegistAtom(Atom & tAtom, int tN);
        
   // private :
        
        SeriNumber   itsSeriNum; 
        ID           itsID;
        
        Residue   *  itsInitRes;
        Residue   *  itsEndRes;
        
        SeriNumber   itsSheetSeriNum; 
        
        Atom      *  itsCurAtom;  // Registration: Atom name in current strand
        Atom      *  itsPrevAtom; // Registration: Atom name in previous strand
        
    };
    
    class Sheet 
    {
    public :
        
        // default constructor 
        Sheet();
  
        // Copy constructor 
        Sheet(const Sheet & tSheet);
        
        // destructor 
        ~Sheet();
        
        inline SeriNumber getSeriNum()
        {
            return itsSeriNum;
        }
        inline void setSeriNum(SeriNumber tN)
        {
            itsSeriNum = tN;
        }
        
        inline ID getID()
        {
            return itsID;
        }
        inline void setID(ID tID)
        {
            itsID = tID;
        }
        
        inline int getNumOfStrands()
        {
            return itsNumOfStrands;
        }
        inline void setNumOfStrands(int tN)
        {
            itsNumOfStrands = tN;
        }
        
        std::vector<Strand>    allStrands;
        std::vector<int>       senses;
        
   //  private :
        
        SeriNumber     itsSeriNum;
        ID             itsID;
        int            itsNumOfStrands;
        
        
    };
    

    

    
}

#endif	/* SECONDARYSTRUCTURES_H */

