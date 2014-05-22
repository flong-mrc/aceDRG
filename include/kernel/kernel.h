//--------------------------------------------------------------------------------
//
//       File Name :   kernel.h  
//       
//       Purpose   :   Interface  
//                     This head file declares 
//                     1) header files for math operations
//                     2) 
// 
//       Date of last modification:   08/07/2011                 
//                                     
//---------------------------------------------------------------------------------


#ifndef KERNEL_H
#define KERNEL_H

#ifndef GLOBAL_H_
#include "global.h"
#endif 

#ifndef CONSTANTS_H_
#include "constants.h"
#endif


namespace LIBMOL
{
    class Entity 
    {
    public:
        
        /* Default constructor, with all member variable initiated
           to the default values
        */
        Entity ()
        {
        }
        
        /* Copy constructor */
        //Entity ( const Entity& tEntity, int tCopylevel );
        
        // Destructor 
        virtual ~Entity ()
        {
        }
        
    private:
        
        
        
        
    };

    
}
#endif


