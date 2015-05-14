#ifndef PrmGroupH
#define PrmGroupH
// =============================================================================
// PrmFileLib - Parametr file parsing and manipulation library
// -----------------------------------------------------------------------------
// (c) 2008 Petr Kulhanek (kulhanek@enzim.hu)
//          Institute of Enzymology, Karolina ut 29, Budapest H-1113, Hungary
// (c) 2008 Martin Petrek (petrek@chemi.muni.cz)
//          National Centre for Biomolecular Research, Kotlarska 2,
//          Brno CZ-611 37, Czech Republic
// -----------------------------------------------------------------------------
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 51 Franklin Street, Fifth Floor,
//    Boston, MA  02110-1301  USA
// =============================================================================

#include "SmallString.hpp"

//------------------------------------------------------------------------------

class CPrmGroup;
class CPrmSection;

//------------------------------------------------------------------------------

class CPrmGroup {
public:
    // constructor and destructor
    CPrmGroup(const CSmallString& name);
    ~CPrmGroup(void);

    //! count unprocessed lines in the group
    int CountULines(void);

    //! count sections in the group
    int CountGroup(void);

    //! return group name
    const CSmallString& GetName(void);

// section of private data ----------------------------------------------------
private:
    CSmallString    Name;           // group name
    CPrmGroup*      Next;           // next group
    CPrmSection*    FirstSection;   // first section in the group
    bool            Processed;      // was the group accessed somehow

    friend class CPrmFile;
};

//------------------------------------------------------------------------------
#endif
