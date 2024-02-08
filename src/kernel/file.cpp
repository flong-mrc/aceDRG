/*
 *
 * File:   file.cpp
 * Author: flong
 *
 * Created on September 22, 2011, 6:17 PM
 */

#include "file.h"

namespace LIBMOL
{
    File::File()
    {

    }


    File::File(std::string tFName, std::ios::openmode tOpenMode=std::ios::in)
    {
        open(tFName.c_str(), tOpenMode);

        itsState.isOpen = is_open();
    }

}
