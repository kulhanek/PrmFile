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

#include <PrmLine.hpp>
#include <PrmFile.hpp>

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CPrmLine::CPrmLine(void)
{
    Next = NULL;
    Processed = false;
}

//------------------------------------------------------------------------------

CPrmLine::~CPrmLine(void)
{

}

//------------------------------------------------------------------------------

bool CPrmLine::IsProcessed(void)
{
    return(Processed);
}

//------------------------------------------------------------------------------

const CSmallString& CPrmLine::GetText(void)
{
    return(Text);
}

//------------------------------------------------------------------------------

bool CPrmLine::Split(CSmallString& key,CSmallString& value)
{
    key   = NULL;
    value = NULL;

    int lto = Text.GetLength() - 1;

    // key starts always on 0 (due to parsing of file)
    int skey = 0;

    // find end of key
    int ekey = Text.Scan(CPrmFile::white_list);

    if( ekey <= 0 ) return(false);       // line is only key, this is not permitted

    // extract key
    key = Text.GetSubStringFromTo(skey,ekey-1);

    // find beggining of value
    int svalue = Text.Verify(CPrmFile::white_list,ekey,lto);

    // value cannot be only empty string because of parsing
    value = Text.GetSubStringFromTo(svalue,lto);

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

