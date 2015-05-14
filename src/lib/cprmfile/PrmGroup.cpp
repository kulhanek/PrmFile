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

#include <PrmGroup.hpp>
#include <PrmSection.hpp>
#include <PrmLine.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CPrmGroup::CPrmGroup(const CSmallString& name)
{
    Name = name;
    Next = NULL;
    FirstSection = NULL;
    Processed = false;
}

//------------------------------------------------------------------------------

CPrmGroup::~CPrmGroup(void)
{
    CPrmSection* p_sec = FirstSection;

    while( p_sec != NULL ) {
        CPrmSection *p_next = p_sec->Next;
        delete p_sec;
        p_sec = p_next;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CPrmGroup::CountULines(void)
{
    CPrmSection*   p_sec = FirstSection;
    int            cnt=0;

    if( Processed == false ) cnt++;

    while( p_sec != NULL ) {
        CPrmLine* p_line = p_sec->FirstLine;
        while( p_line != NULL ) {
            if( p_line->IsProcessed() == false ) cnt++;
            p_line = p_line->Next;
        }
        p_sec = p_sec->Next;
    }

    return(cnt);
}

//------------------------------------------------------------------------------

int CPrmGroup::CountGroup(void)
{
    CPrmSection*   p_sec = FirstSection;
    int            cnt=0;

    while( p_sec != NULL ) {
        cnt++;
        p_sec = p_sec->Next;
    }

    return(cnt);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

const CSmallString& CPrmGroup::GetName(void)
{
    return(Name);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

