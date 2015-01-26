/* 
 * File:   global.h
 * Author: flong
 *
 * Created on August 2, 2011, 1:35 PM
 */

#ifndef GLOBAL_H
#define	GLOBAL_H


//  C++ Standard libs

#ifndef __IOSTREAM_
#include <iostream>
#endif

#ifndef __FSTREAM__
#include <fstream>
#endif

#ifndef __IOMANIP__
#include <iomanip>
#endif

#ifndef  __STRING_
#include <string>
#endif

#ifndef  __CSTRING_
#include <cstring>
#endif

#ifndef __SSTREAM_
#include <sstream>
#endif

#ifndef  __CCTYPE__
#include <cctype>
#endif

#ifndef __CSTDLIB_
#include <cstdlib>
#endif

#ifndef __CMATH_
#include <cmath>
#endif

#ifndef  __CTIME_
#include <ctime>
#endif

//  C standard lib
#ifndef __UNISTD_
#include <unistd.h>
#endif

#ifndef __SYS_STAT_
#include <sys/stat.h>
#endif

//  C++ Standard Template Lib

#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <iterator>
#include <utility>
#include <locale>

// #include <sqlite3.h>

#ifndef DATATYPE_H
#include "datatype.h"
#endif 

#ifndef UTILITY_H
#include "utility.h"
#endif


#endif	/* GLOBAL_H */

