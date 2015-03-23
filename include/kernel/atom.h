/* 
 * File:   atom.h
 * Author: flong
 *
 * Created on August 2, 2011, 6:38 PM
 */

#ifndef ATOM_H
#define	ATOM_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef RESTRAINTLISTS_H
#include "restraintLists.h"
#endif

#ifndef RING_H
#include "ring.h"
#endif

namespace LIBMOL // temp 
{
    class NeighbList;
    class TreeSeg;
    
    class Bond;
    class BondDict;
    class BondList;
    
    class Angle;
    class AngleDict;
    class AngleList;
    
    class Torsion;
    class TorsionDict;
    class TorsionList;
    
    class Chiral;
    class ChiralList;
    
    class RestraintLists;
    
    class Ring;
    class RingDict;
    
/*
 *  Atom class describes the chemical entities, atoms etc. It derived publicly 
 * from class Entity.
 *  @Para:
 *         (1) general properties
 *         name;
 *         elementName;
 *         serNum-- the serial number in an assembly of atoms;
 *         
 *         (2) physical properties 
 *         type  -- a parameter linked to bonding condition, used in 
 *                    building molecules and forcefield optimization;
 *         mass;
 *         coordinates;
 *         coordinate_exp--dictionary values of relative coordinates
 *         coordinates_pre--coordinates in the previous moment;
 *         force;
 *         velocity;
 *         acceleration;  
 *         
 *         (3) PDB and x-ray related
 *          
 *         
 *  
 *         
 * 
 */
   
    enum hasType {
        UNDEFINED,
        DEFINED,
        NOTDETERMINED,
        DETERMINED
    };
    
    class Atom
    {
    public:
        // default  constructor
        Atom();
        // copy constructor
        Atom(const Atom & tAtom);
        
        // destructor 
        ~Atom();
        
        std::string  getName() const;
        void         setName(std::string);
        
        int          getSeriNum() const;
        void         setSeriNum(int tN);
        
        Element      getElementType() const;
        void         setElementType(Element tElem);
        
        IDCode       getID() const;
        void         setIDCode(IDCode tID);
        
        ResName      getResName() const;
        void         setResName(ResName tName);
        
        SeriNumber   getSeqNum() const;
        void         setSeqNum(SeriNumber tN);
        
        InsCode      getInsCode() const;
        void         setInsCode(InsCode tID);
        
        SeriNumber   getSegNum() const;
        void         setSegNum(SeriNumber tN);
        
        IDCode       getSegID() const;
        void         setSegID(IDCode tID);
        
        ChainID      getChainID() const; 
        void         setChainID(ChainID tID);
        
        SeriNumber   getModSeriNum() const;
        void         setModSeriNum(SeriNumber tN);
        
        REAL         getMass() const;
        void         setMass(REAL tMass);
        
        REAL         getCharge() const;
        void         setCharge(REAL tC);
        
        REAL         getPartialCharge() const;
        void         setPartialCharge(REAL tC);
        
        AltLoc       getAltLoc() const;
        void         setAltLoc(AltLoc tA);
        
        REAL         getOccup() const;
        void         setOccup(REAL tO);
        
        REAL         getTempFact() const;
        void         setTempFact(REAL tT);     
        
        REAL         getRadius() const;
        void         setRadius(REAL tR);
        
        REAL         getIonRadius() const;
        void         setIonRadius (REAL tIR);
        
        
        // some data variables made public for efficiency 
        //Cartesian coordinates of a atom in angstroms
        std::vector<REAL> coords; 
        // standard deviation of coordinates
        std::vector<REAL> sigCoords; 
        // the expected values of coordinates
        std::vector<REAL> coords_exp;  
        
        // the Cartesian velocities of the atom 
        std::vector<REAL> velos;   
        // the Cartesian accelerations of the atom
        std::vector<REAL> accels;   
        // the forces acting on the atom
        std::vector<REAL> forces;     
        // sometime atom needs to remember the following  
        std::vector<REAL> coordsPre; 
        std::vector<REAL> velosPre;   
        std::vector<REAL> accelsPre;  
        
        std::vector<REAL> Uxx;      // anisotropic temperature factors
                                    // U[0] = u11, U[1]=u22, U[2]=u33,
                                    // U[3] = u12, U[4]=u13, U[5]=u23  
        std::vector<REAL> sigUxx; 
        
        
        std::vector<REAL>  altLoc;
        std::list<Atom>    neighbAtoms;   
        std::vector<int>   connectedAtoms;
        RestraintLists   * allRestrLists;

        std::list<int>     treeSeg;       // atom serial numbers only
        
        /* possible key-type for map types are 
        *  hydrogen bonding type
        *  ligand type 
        */      
        
        std::map<std::string, int> types;
        
        std::map<std::string, int> state;        
        

    private :
      
        std::string   itsName;
        int           itsSeriNum;
        Element       itsElementType; // the element type of an atom
       
        IDCode        itsID;          // atom role name. See PDB convention
        ResName       itsResName;     // amino acid abbreviation for an atom
        SeriNumber    itsSeqNum;      // Residue sequence number 
        InsCode       itsInsCode;     //  Code for insertion of residues
        SeriNumber    itsSegNum;      // segment serial number
        IDCode        itsSegID;       // segment identifier 
        ChainID       itsChainID;     // the chain ID for an atom 
        
        SeriNumber    itsModNum;       // the serial number of the model
        
        REAL          itsMass;        // the mass of an atom
        REAL          itsCharge;      // the charges linked with atomic number 
        REAL          itsPCharge;     // the partial charge for an atom
 
        REAL          itsRadius;      // atom radius
        REAL          itsIonRadius;   // ionic radius for the atom  
        
        // atoms in macromolecules from x-ray ( e.g. from PDB) 
        AltLoc        itsAltLoc;      // alternate location indicator
        REAL          itsOccup;       // Occupancy for an atom
        
        REAL          itsTempFact;    // B value or temperature factor

        
        // standard deviations for some member variables
        REAL      itsSigOcc;         // standard deviation of occupancy
        REAL      itsSigTemp;        // standard deviation of temperature factor
 
               
        
    };
    
    
    class AtomDict 
    {
    public:
        
        //Default constructor
        AtomDict();
        // Copy constructor
        AtomDict(const AtomDict & tAtom);
        // Destructor 
        ~AtomDict();
        
        int  getNumAtomsWithin2stNB(std::vector<AtomDict> & tAllAtoms);
        int  getNum1stNbHave2edNb(std::vector<AtomDict> & tAllAtoms);
        void setCodClass();
        void outRingSec();
        int  getMinRing();
        int  getMinRing2();
        void setBaseRingProps();
        
        void setCodChemType();
        
        void fromCodClassToEnerType();
              
        void outNeighBAtoms();
        
        int  atomPosition(std::vector<AtomDict> & tAtoms);
        
        // These are from Dictionary cif file
        int      seriNum;
        ID       resName;
        ID       id; 
        ID       chemType;
        Name     enerType;
        REAL     charge;
        REAL     parCharge;
        REAL     formalCharge;
        REAL     radius;
        REAL     ionRadius;
        int      inChiralIdx;
        int      chiralIdx;
        int      bondingIdx;
        bool     isMetal;
        ID       metalGeo;
        REAL     isoB;
        REAL     ocp;
        int      symmMult;
        
        
        // These are CCP4 related 
        ID       ccp4Type;
        
        // These are COD related 
        ID       cChemType;              // COD element ID with sp2 info 
        ID       codClass;
        int      codMolIdx;
        ID       codCifName;
        ID       codAtmRoot;
        ID       codNBSymb;
        ID       codNB2Symb;
        ID       codNB3Symb;
        ID       codAtmMain;
        
       
        int      hashingValue;
        bool     coordExist;
        int      numPi;
        
        bool     isCChemTypeSet;
        bool     isInPreCell;
        bool     chiralChecked;
        
        // symmetry-related 
        ID       sId;
        ID       symmOp;
        // int      bondingIdx;
        /* meaning of curBondingIdx
         * 0      undetermined  
         * 1      sp1
         * 2      sp2
         * 3      sp3
         * 4      sp3d
         * 
         * int    chiralIdx;
         * -1     a plane center
         * 0      undetermined
         * 1      chiral center
         */
        
        
        std::vector<ID> codClassV;
        
        std::map<std::string, int>  existProps;
        // existProps may have the following keys
        // 1. Coords
        // 2. ChemType
        // 3. EnerType;
        // 4. ParCharge;
        
        std::vector<REAL>          coords;
        std::vector<REAL>          fracCoords;
        std::vector<REAL>          forces;
        std::vector<int>           connAtoms;
        std::vector<int>           connHAtoms; // H atoms connected to this atom
        std::vector<int>           neighbAtoms;
        std::vector<int>           inBonds;
        std::vector<AngleDict>     inAngles;
        std::vector<int>           inRings; 
        std::vector<int>           inChirals;
        std::vector< std::string > nbRep;
        std::map<std::string, int> ringRep;
        std::map<std::string, std::string>  ringRepS;
        std::map<std::string, int> ringRepBySeriNum;
        std::map<std::string, std::string> baseRingProp;
        std::map<std::string, std::vector<int> > tree;
        
        REAL                       treeBond;
        REAL                       treeAngle;
        REAL                       treeTorsion;
        
                                            // element 1: immediate neighbor 
                                            // element 2: include neighbor of 
                                            //            of neighbor atoms 
                                            //            => codClass
                                            // element 3: one more neighbor layers
                                            //        than key 2
    };
    
    
    extern int getAtom(std::string             tId,
                       int                     tSeri,
                       std::vector<AtomDict> & tAtoms);
    
    
 
}
#endif	/* ATOM_H */

