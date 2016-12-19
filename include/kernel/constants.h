/* 
 * File:   constants.h
 * Author: flong
 *
 * Created on August 2, 2011, 2:05 PM
 */

#ifndef CONSTANTS_H
#define	CONSTANTS_H

// math constants 

#ifndef __PI__
#define   PI         4.0*atan(1.0)
#endif

#ifndef __PI180__
#define   PI180       PI/180.0
#endif

#ifndef __PID180__
#define   PID180      180.0/3.141592653589 
#endif

#ifndef __GOLD__
#define   GOLD        0.618034        // The gold ratios
#endif

static const int DIM3D = 3; // may not need it 

// constants from chemistry, element types etc



// constants from physics and chemistry

#ifndef __BOLTZ__
#define   BOLTZ       0.86161442E-4  // (1.38047/1.6021E4) 
#endif

// Sone default constants
#define VDWCONST       0.000    // use carbon-carbon value
#define RTHRESHOLD     0.10     // used to check small molecule cif 

// Process related constants

static const int MaxStateStored = 1000;  // number of loacal minimus stored

static const int MaxOptimSets = 100;     

#endif	/* CONSTANTS_H */

