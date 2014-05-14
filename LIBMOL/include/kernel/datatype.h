/*
 * datatype.h
 *
 *  Created on: Jul 11, 2011
 *  last updated: Oct 26, 2011
 *  Author: flong
 */

#ifndef DATATYPE_H
#define DATATYPE_H

// Basic data type defined during "make/cmake" process

#ifndef CONFIG_H
#include "config.h"
#endif 

#ifndef GLOBAL_H
#include "global.h"
#endif 

namespace LIBMOL
{
// numerical types defined for calculations
typedef double REAL;  // selected from float, double and long double 

// general usage
typedef std::size_t Size;
typedef std::string Name;
typedef std::string ID;
typedef std::string ALine;

/* Data types for specific uses  */

typedef  std::string  Date;  // date DD-MMM-YYYY
static const char spac3[5] ="  ";

// for files. paths and entity ID 

typedef const char  *  FileName;
typedef const char  *  SqliteStatment;
typedef char  *  PathName;


// for atom, chain, model and entities defined in a PDB file etc.

typedef int          SeriNumber;   // serial number of an entity 
typedef std::string  IDCode;       // ID code of the entity 
typedef std::string  RecHead;      // Starting word of a record in a file

typedef  std::string SeqDBName;    // sequence database name
typedef  std::string SeqDBAcCode;  // sequence database accession code
typedef  std::string SeqDBIdCode;  // sequence database identification code
typedef  std::string InsCode;      // insertion code
typedef  std::string AltLoc;       // alternate location indicator


typedef  std::string AtomName;     // name of the atom
typedef  std::string Element;      // chemical element name
typedef  std::map<std::string, int> Atomtype; // composite property for a atom

typedef  std::string ResName;      // residue name
typedef  std::vector<std::string>  ListResName; // a vector of residue names

typedef  std::string ChainID;      // ID for a chain 

typedef  std::string HelixID;      // helix ID
typedef  std::string StrandID;     // strand ID
typedef  std::string SheetID;      // sheet ID
typedef  std::string urnID;        // turn ID

// define some common values (usually used as default values
#define ZeroInt    0
#define ZeroShort  0
#define ZeroReal   0.0
#define NNullStrinullChar   ""
#define NullString ""
#define NullPoint  NULL

// May be better using template here
class sortMap
{
  public:
        std::string key;
        int val;
};

struct sortMap2
{
    std::string key;
    int val;
    int nNB;
};

struct sortMap3
{
    std::string key;
    std::string val;
};

class sortLine
{
public:
    std::string key;
    std::string line;
};

}

#endif /* TYPE_H_ */


