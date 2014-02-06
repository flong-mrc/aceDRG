/* 
 * File:   otherAtomContainer.h
 * Author: flong
 *
 * Created on August 11, 2011, 2:57 PM
 */

#ifndef OTHERATOMCONTAINER_H
#define	OTHERATOMCONTAINER_H

#ifndef KERNEL_H
#include "kernel.h"
#endif

namespace LIBMOL
{
    class Atom;
    
    // All classes defined here are simple AtomAssembly or atomContainer
    class Link  : public AtomAssembly
    {
    public :
        
        Link();
        Link(Link & tL);
        virtual ~Link();
        
                
    };
    
    class Helice  : public AtomAssembly
    {
    public :
        
        Helice();
        Helice(Helice & tH);
        virtual ~Helice();
        
                
    };

    class Sheet  : public AtomAssembly
    {
    public :
        
        Sheet();
        Sheet(Sheet & tS);
        virtual ~Sheet();
                
    };
    
    class Turn  : public AtomAssembly
    {
    public :
        
        Turn();
        Turn(Link & tL);
        virtual ~Turn();                       
    };
    
#endif	/* OTHERATOMCONTAINER_H */

