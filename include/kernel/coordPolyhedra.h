/* 
 * File:   coordPolyhedra.h
 * Author: flong
 *
 * Created on October 17, 2012, 10:32 AM
 */

#ifndef COORDPOLYHEDRA_H
#define	COORDPOLYHEDRA_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

#ifndef FILE_H
#include "file.h"
#endif

#ifndef ATOM_H
#include "atom.h"
#endif

#ifndef BOND_H
#include "bond.h"
#endif

#ifndef ANGLE_H
#include "angle.h"
#endif

#ifndef TORSION_H
#include "torsion.h"
#endif

#ifndef RESIDUE_H
#include "residue.h"
#endif

#ifndef RING_H
#include "ring.h"
#endif

#ifndef PLANE_H
#include "plane.h"
#endif

#ifndef CHAIN_H
#include "chain.h"
#endif

#ifndef SSBOND_H
#include "ssbond.h"
#endif

#ifndef LIBG_LINK_H
#include "libglink.h"
#endif

#ifndef UTILITY_H
#include "utility.h"
#endif

#ifndef PERIODICTABLE_H
#include "periodicTable.h"
#endif


namespace LIBMOL
{
    class Atom;
    class AtomDict;
    
    class Bond;
    class BondDict;
    
    class Ring;
    class RingDict;
    
    class Angle;
    class AngleDict;
    
    class Torsion;
    class TorsionDict;
    
    class Chiral;
    class ChiralDict;
    
    class Plane;
    class PlaneDict;
    
    class PeriodicTable;
    
    class coordPolyhedra 
    {
    public :
        
        // Default constructor
        coordPolyhedra();
        
        // Copy constructor
        coordPolyhedra(coordPolyhedra & tPh);
        
        // Constructor using a coordination number, a group of atoms,
        // and a set of bonds
        coordPolyhedra(int tCoordNum, 
                       std::vector<AtomDict> & tAtoms,
                       std::vector<BondDict> & tBonds);
        
        // Destructor 
        ~coordPolyhedra();
        
        // select a polyhedra and set the default geometry (atom coordinates) 
        // according to the coordination number, atom types and bonds
        void setupGeometry();
        
        // coordination number 2
        void setupLine();         
        
        // coordination number 3
        void setupTriPlanar();    
        void setupTriPyramid();
        void setupTShape();
        
        // coordination number 4
        void setupTetraHedral();
        void setupSqPlanar();
        
        // coordination number 5
        void setupSqPyramid();
        void setupTriBipyramid();
        
        // coordination number 6
        void setupHexPlanar();
        void setupTriPrism();
        void setupOctahedral();                // O_h
        
        // coordination number 7
        void setupCappedOctahedron();          // C_3v
        void setupCappedtriPrism();            // C_2v
        void setuoPentaBipyramid();            // D_5h
        
        
        // coordination number 8     
        void setupDodecahedron();              // D_2d
        void setupCube();                      // O_h
        void setupSqAntiprism();               // D_4d 
        void setupHexBipyramid();              // D_6h
        
        // coordination number 9
        void setup3FaceCenTriPrism();          // D_3h
        
        // coordination number 10
        void setupBicappedSqAntiprism();       // D_4d
        
        // coordination number 11 
        void setupAllFacedCappedTriPrism();    // D_3h
        
        // Coordination Number 12
        void setupCubeOctahedron();            // O_h
        
        
        // permute ligand atoms within a polyhedra 
        
        
        int                         coordNum;
        
        // Convention: 
        // the first atom is metal one if we deal with a metal cluster
        std::vector<AtomDict>       atoms;   
        
        std::vector<BondDict>       bonds;
        
        
    };
    
    /* The following class generating geometrical parameters (bonds, angles) 
     * based on input coordination numbers and geometric key words.
     * 
     *  
     */
    class CoordGeom 
    {
    public:
        
        // Default
        CoordGeom();
        
        // Copy constructor
        CoordGeom(CoordGeom & aCGO);
        
        // Destructor 
        ~CoordGeom();
       
        void genBondAndAngles();
        
        void genBondAndAnglesLINEAR();
        
        
        int             coordNum;
        std::string     coordGeom;
        
        std::vector<BondDict>               bondCandidates;
        std::vector<AngleDict>              angleCandidates;
       
        
        
        
        
    };
    
}    

#endif	/* COORDPOLYHEDRA_H */

