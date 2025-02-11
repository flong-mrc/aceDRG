/* 
 * File:   DnaRna.h
 * Author: flong
 *
 * Created on December 6, 2011, 2:21 PM
 */

#ifndef DNARNA_H
#define	DNARNA_H

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

#ifndef SECONDARYSTRUCTURES_H
#include "secondaryStructures.h"
#endif

namespace LIBMOL
{
    class Atom;
    class Residue;
    class ModRes;
    class Chain;

    class Helix;
    class Sheet;
    class Strand;
    class Turn;
    
    class DnaBase 
    {
    public:
        
        // default constructor
        DnaBase();
        // copy constructor
        DnaBase(const DnaBase & tBase);
        
        // destructor
        ~DnaBase();
        
        
        std::vector<Atom>        atoms;
        
        std::string              baseID;    // one of A, G, T, C, U 
        std::string              fullName;  // one of Adenine, Guanine, 
                                            // Thymine, Cytosine and Uracil
        std::string              group;     // Purine or Pyrimidine 
        std::string              cID;
        
    };
    
    class DnaBasePair 
    {
    public:
        // default constructor
        DnaBasePair();
        // copy constructor
        DnaBasePair(const DnaBasePair & tBPair);
        // destructor 
        ~DnaBasePair();
        
        void                   setXyzAxis();
        REAL                   calcPTwist();
        REAL                   calcBuckle();
        REAL                   calcOpen();
        REAL                   calcShear();
        REAL                   calcStretch();
        REAL                   calcStagger();
        
        
        std::vector<DnaBase>   bases;
        std::vector<Bond>      hydroBonds;
        
        // parameters describe geometry of the base-pair
        std::vector<REAL>      xyzAxis;
        REAL                   pTwist;
        REAL                   buckle;
        REAL                   opening;
        REAL                   shear;
        REAL                   stretch;
        REAL                   stagger;  
    
    };
    
    // Two successive base-pairs
    class Dinucleotide 
    {
    public :
        // Default constructor
        Dinucleotide();
        // copy constructor 
        Dinucleotide(Dinucleotide & tD);
        // destructor 
        ~Dinucleotide();
        
        
        void setDiXyzAxis();
        void setHelXyzAxis();
        
        
        REAL calcShift();
        REAL calcTilt();
        REAL calcSlide();
        REAL calcRoll();
        REAL calcRise();
        REAL calcTwist();
        
        REAL calcDx();
        REAL calcDy();
        REAL calcInclin();
        REAL calcTip();
        
        // two component base-pairs  
        std::vector<DnaBasePair>     basePairs;
        // parameters describe The relative position and 
        // orientation of successive base pairs
        
        // a dimer reference frame
        std::vector<REAL>            dimeXYZAxis;
        REAL                         shift;
        REAL                         tilt;
        REAL                         slide;
        REAL                         roll;
        REAL                         rise;
        REAL                         twist;
        
        // a local helical frame
        std::vector<REAL>            helixXYZAxis;
        REAL                         dX;
        REAL                         dY;
        REAL                         inclin;
        REAL                         tip;
              
        
    };
    
}

#endif	/* DNARNA_H */

