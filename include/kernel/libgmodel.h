/* 
 * File:   libgmodel.h
 * Author: flong
 *
 * Created on August 11, 2011, 2:22 PM
 */

#ifndef LIBG_MODEL_H
#define	LIBG_MODEL_H

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

#ifndef LIBG_LINK_H
#include "libglink.h"
#endif

#ifndef SECONDARYSTRUCTURES_H
#include "secondaryStructures.h"
#endif

namespace LIBMOL // temp 
{
    class Chain;
    class Residue;
    class Atom;
    class Molecule; 
    class Link;
    class Helix;
    class Sheet;
    
    
    // class Model represents a model of protein.
    
    class Model 
    {
    public:
        
        // default constructor
        Model();
        // copy constructor 
        Model(const Model & tM);
        // constructor using mini info
        Model(ID tID, SeriNumber tN, Name tName);
        
        // destructor 
        virtual ~Model();
 
 
        ID getID() const;
        void setID(ID tID);
        SeriNumber getSeriNum() const;
        void setSeriNum(SeriNumber tSer);
        SeriNumber getModSeriNum();
        void  setModSeriNum(SeriNumber tN);
        
        Model * create(ID tID, 
                       SeriNumber tN, 
                       Name tName,
                       std::vector<Chain> tChains);
        
        Size     getNumOfChains() const;  
        Chain &  getOneChain(SeriNumber tI);
        void     addOneChain(const Chain & tC);
        bool     deleteOneChain(SeriNumber tN);
        
        Size     getNumOfLinks() const;  
        Link  &  getOneLink(SeriNumber tN);
        void     addOneLink(const Link & tL);
        bool     deleteOneLink(SeriNumber tN);
        
       
        Size     getNumOfHelices() const;  
        Helix &  getOneHelix(SeriNumber tI);
        void     addOneHelix(const Helix & tH);
        void     deleteOneHelix(SeriNumber tN);
          
        Size     getNumOfSheets() const;
        Sheet &  getOneSheet(SeriNumber tI);
        void     addOneSheet(const Sheet & tS);
        bool     deleteOneSheet(SeriNumber tN);
        
       // Size     getNumOfTurns();  
       // Turn   & getOneTurn(SeriNumber tN);
       //  void     addOneTurn(Turn & tL);
       //  bool     deleteOneTurn(SeriNumber tN);
        
        
        std::map<std::string, bool>     state;
        
        std::vector<Chain>        chains;
        std::vector<Link>         links;  
        std::vector<Helix>          helices;
        std::vector<Sheet>          sheets;
        
    private :
        ID                        itsID;
        SeriNumber                itsSeriNum;
        

        // std::list<Turn>           itsTurns;
        
    };
    
}


#endif	/* MODEL_H */

