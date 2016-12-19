/* 
 * File:   spaceGroupTable.h
 * Author: flong
 *
 * Created on August 7, 2013, 6:46 PM
 */

#ifndef CRYSTINFO_H
#define	CRYSTINFO_H

#ifndef KERNEL_H
#include "kernel.h"
#endif


namespace LIBMOL
{
 
    class SpaceGroupMember
    {
    public:
        // default constructor
        SpaceGroupMember();
        
        // copying constructor
        SpaceGroupMember(const SpaceGroupMember & aSG);
        
        // destructor
        ~SpaceGroupMember();
        
        int                                                        sgNum;
        std::map<ID, std::vector<ID> >                             sgSymb;
        std::map<std::string, std::vector<std::vector<REAL> > >    sgOp;  // 4 x 4 for each op
        
    };
    
    class SpaceGroupTable
    {
        
    public:
        
        // Default constructor
        SpaceGroupTable();
        
        // default destructor
        ~SpaceGroupTable();
        
        
        void getSgtableFromCCP4();
       
        std::vector<SpaceGroupMember>          allSgs;
        
    private:
        
        SpaceGroupMember     *                 itsCurSG;
        
    };
    
    class Cell 
    {
    public:
        //Default constructor
        Cell();
        // copy constructor
        Cell(const Cell & tC);
        
        //Default destructor
        ~Cell();
        
        void cellOrthoFromFract();
        void cellFractFromOrtho();
        
        
        REAL                                a;
        REAL                                b;
        REAL                                c;
        REAL                                alpha;
        REAL                                beta;
        REAL                                gamma;
        
        REAL                                vol;
        
        std::string                         lattice;  
        
    };
    
    class Resolution
    {
    public :
        // Default constructor
        Resolution();
        
        // Copy constructor
        Resolution(const Resolution & tResol);
        
        // Destructor 
        ~Resolution();
        
        void setResol();
        
        REAL                          resolLimit;
        REAL                          dMax;
        REAL                          thetaMax;
        REAL                          wavLen;
        bool                          lSet;
        
        
        
    };
    
     /* Class crystInfo contains variables and function that describe the 
     * crystal symmetrical properties.
     */
    
    class CrystInfo 
    {
    public :
        
        // default constructor 
        CrystInfo();
        
        // copy constructor
        CrystInfo(const CrystInfo & tCrystInfo); 
        
        
        // destructor
        ~CrystInfo();
        
        
        std::vector<REAL>         ORIGXn;
         
        std::vector<REAL>         SCALEn;
        
        std::vector<REAL>         MTRIXn;       
        
        Cell                 *    itsCell;
        SpaceGroupMember     *    itsSpaceGroup;
        Resolution           *    itsResolution;

    };
    
 
    
    
}
#endif	/* SPACEGROUPTABLE_H */

