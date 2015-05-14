#ifndef PrmSectionH
#define PrmSectionH
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

class CPrmSection;
class CPrmLine;

//------------------------------------------------------------------------------

class CPrmSection {
public:
    // constructor and destructor
    CPrmSection(const CSmallString& name);
    ~CPrmSection(void);

    //! count unprocessed lines in the group
    int CountULines(void);

    //! return section name
    const CSmallString& GetName(void);

    //! count number of lines in the section
    int CountSection(void);

// section of private data ----------------------------------------------------
private:
    CSmallString    Name;       // section name
    CPrmSection*    Next;       // next section
    CPrmLine*       FirstLine;  // first line in the section
    bool            Processed;  // was the section accessed somehow

    friend class CPrmFile;
    friend class CPrmGroup;
};

//------------------------------------------------------------------------------
#endif
