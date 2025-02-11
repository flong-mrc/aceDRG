/* 
 * File:   crysInfo.h
 * Author: flong
 *
 * Created on August 31, 2011, 11:03 AM
 */

#ifndef LIBG_CRYSTINFO_H
#define	LIBG_CRYSTINFO_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

namespace LIBMOL
{
    
    struct  Cell 
    {
        
        REAL        a;
        REAL        b;
        REAL        c;
        
        REAL        alpha;
        REAL        beta;
        REAL        gamma;
        
    };
    /* Class crystInfo contains variables and function that describe the 
     * crystal symmetrical properties.
     */
    
    class CrystInfo 
    {
    public :
        
        // default constructor 
        inline CrystInfo() : itsCell(NullPoint),
                itsSpaceGroup(NullString),
                itsZ(ZeroInt)
        {           
        }
        
        // copy constructor
        inline CrystInfo(CrystInfo & tCrystInfo) 
        {
            itsCell         = tCrystInfo.getCell();
            itsSpaceGroup   = tCrystInfo.getSpaceGroup();
            itsZ            = tCrystInfo.getZ();
        }
        
        // Detailed constructor
        inline CrystInfo(REAL tA, REAL tB, REAL tC,
                  REAL tAl, REAL tBe, REAL tGa,
                  std::string tS, int tZ)
        {
            itsCell = new Cell();
            setCell(tA, tB, tC, tAl, tBe, tGa);
            setSpaceGroup(tS);
            setZ(tZ);
            
        }
        
        // destructor
        inline ~CrystInfo()
        {
            if(itsCell)
            {
                delete itsCell;
                itsCell = NullPoint;
            }
        }
        
        inline Cell *  getCell() const
        {
            return itsCell; 
        }
        inline void setCell(REAL tA, REAL tB, REAL tC,
                            REAL tAl, REAL tBe, REAL tGa) 
        {
            itsCell->a = tA;
            itsCell->b = tB;
            itsCell->c = tC;
            
            itsCell->alpha = tAl;
            itsCell->beta  = tBe;
            itsCell->gamma = tGa;
        }
        
        inline REAL getCell_a() const
        {
            return itsCell->a;
        }
        inline void setCell_a(REAL tA)
        {
            itsCell->a = tA;
        }
           
        inline REAL getCell_b() const
        {
            return itsCell->b;
        }         
        inline void setCell_b(REAL tB)
        {
            itsCell->b = tB;
        }
        
        inline REAL getCell_c() const
        {
            return itsCell->c;
        }         
        inline void setCell_c(REAL tC)
        {
            itsCell->c = tC;
        }      
        
        inline REAL getCell_alpha() const
        {
            return itsCell->alpha;
        }         
        inline void setCell_alpha(REAL tA)
        {
            itsCell->alpha = tA;
        } 
        
        inline REAL getCell_beta() const
        {
            return itsCell->beta;
        }         
        inline void setCell_beta(REAL tB)
        {
            itsCell->beta = tB;
        }
        
        inline REAL getCell_gamma() const
        {
            return itsCell->gamma;
        }         
        inline void setCell_gamma(REAL tG)
        {
            itsCell->gamma = tG;
        }
        
        inline std::string getSpaceGroup () const
        {
            return itsSpaceGroup;
        }
        inline void setSpaceGroup(std::string tS)
        {
            itsSpaceGroup = tS;
        }
        
        inline int   getZ() const
        {
            return itsZ;
        }
        inline void setZ(int tZ)
        {
            itsZ  = tZ;
        }
        
    std::vector<REAL>         ORIGXn;
         
    std::vector<REAL>         SCALEn;
        
    std::vector<REAL>         MTRIXn;       
        
    private :
        
        Cell     *                itsCell;
        std::string               itsSpaceGroup;
        int                       itsZ;
        


    };
    
    
}


#endif	/* CRYSINFO_H */

