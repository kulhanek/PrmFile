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

#include <PrmFile.hpp>
#include <PrmGroup.hpp>
#include <PrmSection.hpp>
#include <PrmLine.hpp>
#include <PrmUtils.hpp>
#include <stdlib.h>
#include <ErrorSystem.hpp>
#include <errno.h>
#include <string.h>
#include <TerminalStr.hpp>

//------------------------------------------------------------------------------

using namespace std;

//------------------------------------------------------------------------------

const char* CPrmFile::white_list = "\t \n";
const char* CPrmFile::comm_list = "!*#";
const char* CPrmFile::quot_list = "'\"";

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

CPrmFile::CPrmFile(void)
{
    _FirstGroup = NULL;
    CurrentGroup = NULL;
    CurrentSection = NULL;
    CurrentLine = NULL;
    FieldPosition = 0;
}

//------------------------------------------------------------------------------

CPrmFile::~CPrmFile(void)
{
    Clear();
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPrmFile::Clear(void)
{
    Name = NULL;
    CurrentGroup = NULL;
    CurrentSection = NULL;
    CurrentLine = NULL;
    FieldPosition = 0;

    CPrmGroup* p_grp = _FirstGroup;

    while( p_grp != NULL ) {
        CPrmGroup* p_ngrp = p_grp->Next;
        delete p_grp;
        p_grp = p_ngrp;
    }

    _FirstGroup = NULL;
}

//------------------------------------------------------------------------------

bool CPrmFile::Read(const CSmallString& name)
{
    // clear previous data
    Clear();

    // open file
    FILE* p_file = fopen(name,"r");
    if( p_file == NULL ) {
        CSmallString    error;
        error << "unable to open file '" << name << "' (" << strerror(errno) << ")";
        ES_ERROR(error);
        return(false);
    }

    Name = name;

    bool result = true;

    // parse file
    if( Parse(p_file) == false ) {
        Clear();
        result = false;
    }

    // close file
    fclose(p_file);

    return(result);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CPrmFile::Parse(FILE* fin)
{
    CSmallString sline;
    int          lineid = 0;

    while( sline.ReadLineFromFile(fin) != false ) {
        lineid++;
        // parse line -------------------
        if( ParseLine(lineid,sline) == false ) {
            CSmallString error;
            error << "unable to parse line number " << lineid << " (" << sline << ")";
            ES_ERROR(error);
            return(false);
        }
        sline = NULL;
    }

    // try to open the first group
    OpenGroup("MAIN");

    return(true);
}

//------------------------------------------------------------------------------

bool CPrmFile::ParseLine(int lineid,const CSmallString& sline)
{
    int tstart, tend;

    // skip if empty line
    tstart = sline.Verify(white_list);
    if( tstart == -1 ) return(true);                   // it is blank line

    // determine type of line
    switch (sline[tstart]) {
    case '!':
    case '#':
    case '*':
        // it is comment skip it
        return(true);
    case '{': {
        // it is group
        tend = tstart + sline.IndexOf('}',tstart);
        if( tend <= tstart+1 ) {
            AddError(lineid,"Group name is not terminated or it is empty string!\n");
            return(false);
        }
        CSmallString    gname;
        gname = sline.GetSubString(tstart+1,tend-1);
        return( AddGroup(lineid,gname) );
    }

    case '[': {
        // it is section
        tend = tstart + sline.IndexOf(']',tstart);
        if( tend <= tstart+1 ) {
            AddError(lineid, "Section name is not terminated or it is empty string!\n");
            return(false);
        }
        CSmallString    sname;
        sname = sline.GetSubString(tstart+1,tend-1);
        return( AddSection(lineid,sname) );
    }

    default:
        // it is line
        return( AddLine(lineid,sline) );
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CPrmFile::AddGroup(int lineid,const CSmallString& gname)
{
    CPrmGroup* p_grp;

    try {
        p_grp = new CPrmGroup(gname);
    } catch(...) {
        return(false);
    }

    if( _FirstGroup == NULL ) {
        _FirstGroup = p_grp;
        CurrentGroup = _FirstGroup;
    } else {
        if( CurrentGroup == NULL ) {
            delete p_grp;
            return(false);
        }
        CurrentGroup->Next = p_grp;
        CurrentGroup = CurrentGroup->Next;
    }

    CurrentSection = NULL;
    CurrentLine = NULL;
    FieldPosition = 0;

    return(true);
}

//------------------------------------------------------------------------------

bool CPrmFile::AddSection(int lineid,const CSmallString& sname)
{
    if( CurrentGroup == NULL ) {
        if( AddGroup(lineid,PRMFILE_MAIN_GROUP_NAME) == false ) return(false);
    }

    CPrmSection* p_sec;

    try {
        p_sec = new CPrmSection(sname);
    } catch(...) {
        return(false);
    }

    if( CurrentGroup->FirstSection == NULL ) {
        CurrentGroup->FirstSection = p_sec;
        CurrentSection = CurrentGroup->FirstSection;
    } else {
        CurrentSection->Next = p_sec;
        CurrentSection = CurrentSection->Next;
    }

    CurrentLine = NULL;
    FieldPosition = 0;

    return(true);
}

//------------------------------------------------------------------------------

bool CPrmFile::AddLine(int lineid,const CSmallString& line)
{
    // beggining
    int tstart = line.Verify(white_list);

    // try to detect beginning of comment
    int pos = tstart;
    int tend = tstart;
    int llen = line.GetLength()-1;

    char fquot[2];
    fquot[0]='\0';
    fquot[1]='\0';

    while( pos <=  llen) {
        // skip quotation
        int quot = line.Scan(quot_list,pos,llen);
        if( quot != -1 ){
            fquot[0] = line[quot];
            pos = quot + 1;
            if( pos > llen ){
                tend = llen;
                break;
            }
            quot = line.Scan(fquot,pos,llen);
            if( quot != -1 ){
                tend = quot;
                break;
            } else {
                tend = llen;
                break;
            }
        }

        // skip comment
        pos = line.Scan(comm_list,pos,llen);
        if( pos == -1 ) {   // if no comment or comment in the current column
            tend = line.GetLength() - 1;
            break;
        }
        tend = pos;
        // check if it is real comment
        // if previous character is not white space then it is not comment
        if( line.Scan(white_list,tend-1,tend-1) == tend-1 ) {
            tend = tend - 1;
            break;
        }
        tend = tend + 1;
        pos = tend;
    }

    // now remove spaces between line end and the remaining part
    tend = line.Verify(white_list,0,tend,true);

    // check if there is active section
    if( CurrentSection == NULL ) {
        AddError(lineid,">>> ERROR: None section is active to add this line!");
        return(false);
    }

    // allocate new line
    CPrmLine* p_line;

    try {
        p_line = new CPrmLine;
    } catch(...) {
        return(false);
    }

    // init section
    p_line->Text = line.GetSubStringFromTo(tstart,tend);
    p_line->Processed = false;
    p_line->Next = NULL;

    // and connect it to current section
    if( CurrentLine != NULL ) {
        CurrentLine->Next = p_line;
        CurrentLine = CurrentLine->Next;
    } else {
        CurrentSection->FirstLine = p_line;
        CurrentLine = CurrentSection->FirstLine;
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

void CPrmFile::AddError(int lineid,const CSmallString& err)
{
    fprintf(stderr,">>> ERROR: %s",(const char*)err);
    if( lineid > 0 ) {
        fprintf(stderr,"           Name : %s\n",(const char*)Name);
        fprintf(stderr,"           Line : %d\n",lineid);
    }
}


//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CPrmFile::Dump(FILE* fout,bool unprocessed)
{
    CTerminalStr str;
    if( fout == NULL ){
        str.Attach(stdout);   // use std output
    } else {
        str.Attach(fout);
    }

    return( Dump(str,unprocessed) );
}

//------------------------------------------------------------------------------

bool CPrmFile::DumpGroup(const CSmallString& gname,FILE* fout,bool unprocessed)
{
    CTerminalStr str;
    if( fout == NULL ){
        str.Attach(stdout);   // use std output
    } else {
        str.Attach(fout);
    }

    return( DumpGroup(gname,str,unprocessed) );
}

//------------------------------------------------------------------------------

bool CPrmFile::Dump(ostream& out,bool unprocessed)
{
    CPrmGroup* p_grp = _FirstGroup;

    while( p_grp != NULL ) {
        if( DumpGroup(p_grp,out,unprocessed) == false ) return(false);
        p_grp = p_grp->Next;
    }

    return(true);
}

//------------------------------------------------------------------------------

bool CPrmFile::DumpGroup(const CSmallString& gname,ostream& out,bool unprocessed)
{
    CPrmGroup*     p_grp = _FirstGroup;
    CSmallString   ugname,xgname;

    // find group with name
    ugname = gname;
    ugname.ToUpperCase();

    while( p_grp != NULL ) {
        xgname = p_grp->GetName();
        xgname.ToUpperCase();
        if( xgname == ugname ) break;
        p_grp = p_grp->Next;
    }

    if( p_grp != NULL ) return( DumpGroup(p_grp,out,unprocessed) );

    return(false);
}

//------------------------------------------------------------------------------

bool CPrmFile::DumpGroup(CPrmGroup* p_grp,ostream& out,bool unprocessed)
{
    if( p_grp == NULL ) return(false);

    if( (unprocessed == true) && (p_grp->CountULines() == 0) ) return(true);  // nothing to print

    // print group name
    out << "{" << p_grp->GetName() << "}" << endl;

    // loop over sections
    CPrmSection* p_sec = p_grp->FirstSection;

    while( p_sec != NULL ) {
        if( (unprocessed == true) && (p_sec->CountULines() == 0) ) {
            p_sec = p_sec->Next;
            continue;   // nothing to print in this section
        }

        // print section name
        out << "    [" << p_sec->GetName() << "]" << endl;

        // loop over lines
        CPrmLine* p_line = p_sec->FirstLine;

        while( p_line != NULL ) {
            if( unprocessed == true ) {
                if( p_line->IsProcessed() == false ){
                    out << "        " << p_line->GetText() << endl;
                }
            } else {
                out << "        " << p_line->GetText() << endl;
            }
            p_line = p_line->Next;
        }

        p_sec = p_sec->Next;
    }

    return(true);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

int CPrmFile::CountULines(void)
{
    CPrmGroup* p_group = _FirstGroup;
    int        counter = 0;

    while( p_group != NULL) {
        counter += p_group->CountULines();
        p_group = p_group->Next;
    }

    return(counter);
}

//------------------------------------------------------------------------------

int CPrmFile::CountULines(const CSmallString& gname)
{
    CSmallString   ugname,xgname;

    // find group ----------------------------------
    ugname = gname;
    ugname.ToUpperCase();

    CPrmGroup *p_group = _FirstGroup;
    while( p_group != NULL) {
        xgname = p_group->GetName();
        xgname.ToUpperCase();
        if( ugname == xgname ) break;
        p_group = p_group->Next;
    }

    // count line or return -1 if group was not found
    if( p_group != NULL ) return( p_group->CountULines() );

    return(-1);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CPrmFile::FirstGroup(void)
{
    CurrentGroup = _FirstGroup;
    CurrentSection = NULL;
    CurrentLine = NULL;
    FieldPosition = 0;

    if( CurrentGroup ){
        CurrentGroup->Processed = true;
    }

    return( CurrentGroup != NULL );
}

//------------------------------------------------------------------------------

bool CPrmFile::OpenGroup(const CSmallString& gname)
{
    if( ! FirstGroup() ) return(false);

    CSmallString ugname,xgname;

    ugname = gname;
    ugname.ToUpperCase();

    while( CurrentGroup != NULL ) {
        xgname = CurrentGroup->GetName();
        xgname.ToUpperCase();
        if( xgname == ugname ){
            CurrentGroup->Processed = true;
            return(true);
        }
        CurrentGroup = CurrentGroup->Next;
    }

    return(false);
}

//------------------------------------------------------------------------------

bool CPrmFile::NextGroup(void)
{
    if( CurrentGroup != NULL ) CurrentGroup = CurrentGroup->Next;

    CurrentSection = NULL;
    CurrentLine = NULL;
    FieldPosition = 0;

    if( CurrentGroup != NULL ) {
        CurrentSection = CurrentGroup->FirstSection;
        if( CurrentSection != NULL ) {
            CurrentLine = CurrentSection->FirstLine;
        }
    }

    if( CurrentGroup ){
        CurrentGroup->Processed = true;
    }

    return( CurrentGroup != NULL );
}

//------------------------------------------------------------------------------

int CPrmFile::CountGroup(void)
{
    if( CurrentGroup != NULL ) return( CurrentGroup->CountGroup() );
    return(-1);
}

//------------------------------------------------------------------------------

const CSmallString& CPrmFile::GetGroupName(void)
{
    static CSmallString null;

    if( CurrentGroup != NULL ) {
        return( CurrentGroup->GetName() );
    } else {
        return( null );
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CPrmFile::FirstSection(void)
{
    CurrentSection = NULL;
    CurrentLine = NULL;
    FieldPosition = 0;

    if( CurrentGroup != NULL ) CurrentSection = CurrentGroup->FirstSection;

    if( CurrentSection ){
        CurrentSection->Processed = true;
    }

    return( CurrentSection != NULL );
}

//------------------------------------------------------------------------------

bool CPrmFile::OpenSection(const CSmallString& sname)
{
    if( ! FirstSection() ) return(false);

    CSmallString usname,xsname;

    usname = sname;
    usname.ToUpperCase();

    while( CurrentSection != NULL ) {
        xsname = CurrentSection->GetName();
        xsname.ToUpperCase();
        if( xsname == usname ){
            CurrentSection->Processed = true;
            return(true);
        }
        CurrentSection = CurrentSection->Next;
    }

    return(false);
}

//------------------------------------------------------------------------------

int CPrmFile::OpenSectionAndCount(const CSmallString& sname)
{
    if( ! FirstSection() ) return(-1);

    CSmallString usname,xsname;

    usname = sname;
    usname.ToUpperCase();

    while( CurrentSection != NULL ) {
        xsname = CurrentSection->GetName();
        xsname.ToUpperCase();
        if( xsname == usname ) {
            CurrentSection->Processed = true;
            return(CurrentSection->CountSection());
        }
        CurrentSection = CurrentSection->Next;
    }

    return(-1);
}

//------------------------------------------------------------------------------

bool CPrmFile::NextSection(void)
{
    if( CurrentSection != NULL ) CurrentSection = CurrentSection->Next;

    CurrentLine = NULL;
    FieldPosition = 0;

    if( CurrentSection != NULL ) CurrentLine = CurrentSection->FirstLine;

    if( CurrentSection ){
        CurrentSection->Processed = true;
    }

    return( CurrentSection != NULL );
}

//------------------------------------------------------------------------------

int CPrmFile::CountSection(void)
{
    if( CurrentSection != NULL ) return(CurrentSection->CountSection());
    return(-1);
}

//------------------------------------------------------------------------------

const CSmallString& CPrmFile::GetSectionName(void)
{
    static CSmallString null;

    if( CurrentSection != NULL ) {
        return( CurrentSection->GetName() );
    } else {
        return( null );
    }
}

//------------------------------------------------------------------------------

void CPrmFile::SetSecAsProcessed(void)
{
    CPrmLine* p_line = NULL;
    if( CurrentSection != NULL ) p_line = CurrentSection->FirstLine;

    if( CurrentSection ){
        CurrentSection->Processed = true;
    }

    while( p_line != NULL ) {
        p_line->Processed = true;
        p_line = p_line->Next;
    }
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CPrmFile::FirstLine(void)
{
    CurrentLine = NULL;
    FieldPosition = 0;

    if( CurrentSection != NULL ) CurrentLine = CurrentSection->FirstLine;
    return( CurrentLine != NULL );
}

//------------------------------------------------------------------------------

bool CPrmFile::NextLine(void)
{
    FieldPosition = 0;

    if( CurrentLine != NULL ) CurrentLine = CurrentLine->Next;
    return( CurrentLine != NULL );
}

//------------------------------------------------------------------------------

bool CPrmFile::GetLineKey(CSmallString& key)
{
    if( CurrentLine == NULL ) return(false);

    CSmallString value;
    return( CurrentLine->Split(key,value) );
}

//------------------------------------------------------------------------------

bool CPrmFile::GetLineValue(CSmallString& value)
{
    if( CurrentLine == NULL ) return(false);

    CSmallString key;
    return( CurrentLine->Split(key,value) );
}

//------------------------------------------------------------------------------

bool CPrmFile::GetLine(CSmallString& value)
{
    if( CurrentLine == NULL ) return(false);
    value = CurrentLine->Text;
    return(true);
}

//------------------------------------------------------------------------------

void CPrmFile::SetCurrentLineProcessed(void)
{
    if( CurrentLine != NULL ) CurrentLine->Processed = true;
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CPrmFile::GetIntegerByKey(const CSmallString& name, int& value)
{
    CSmallString svalue;
    int          lvalue;
    char*        endptr = NULL;

    // find key --------------------------------------
    if( FindKey(name,svalue) == true ) {
        if( svalue == NULL ) return(false);
        // try to decode data
        lvalue=strtol(svalue, &endptr, 10);
        if( endptr == NULL ) return(false);
        if( *endptr == '\0' ) {
            value = lvalue;
            SetCurrentLineProcessed();
            return(true);
        }
    }

    return(false);
}

//------------------------------------------------------------------------------

bool CPrmFile::GetStringByKey(const std::string& name, std::string& value)
{
    CSmallString iname = name;
    CSmallString ivalue = value;
    bool ret = GetStringByKey(iname,ivalue);
    value = ivalue;
    return(ret);
}

//------------------------------------------------------------------------------

bool CPrmFile::GetStringByKey(const CSmallString& name, CSmallString& value)
{
    CSmallString svalue;

    // find key --------------------------------------
    if( FindKey(name,svalue) == false )  return(false);

    // get string ------------------------------------
    value = svalue;

    // strip down quotation marks if they are present
    int sbeg = 0;
    int send = svalue.GetLength()-1;
    if( send > 1 ) {
        if( (value[0] == '\'') && (value[send] == '\'')  ) {
            sbeg = sbeg + 1;
            send = send - 1;
        } else if( (value[0] == '"') && (value[send] == '"') ) {
            sbeg = sbeg + 1;
            send = send - 1;
        }
    }

    if( sbeg <= send ) {
        value = svalue.GetSubStringFromTo(sbeg,send);
    } else {
        value = "";
    }

    SetCurrentLineProcessed();
    return(true);
}

//------------------------------------------------------------------------------

bool CPrmFile::GetFloatByKey(const CSmallString& name, float& value)
{
    CSmallString svalue;
    float        lvalue;
    char*        endptr = NULL;

    // find key --------------------------------------
    if( FindKey(name,svalue) == true ) {
        if( svalue == NULL ) return(false);
        // try to decode data
        lvalue=strtof(svalue, &endptr);
        if( endptr == NULL ) return(false);
        if( *endptr == '\0' ) {
            value = lvalue;
            SetCurrentLineProcessed();
            return(true);
        }
    }

    return(false);
}

//------------------------------------------------------------------------------

bool CPrmFile::GetDoubleByKey(const CSmallString& name, double& value)
{
    CSmallString svalue;
    double       lvalue;
    char*        endptr = NULL;

    // find key --------------------------------------
    if( FindKey(name,svalue) == true ) {
        if( svalue == NULL ) return(false);
        // try to decode data
        lvalue=strtod(svalue, &endptr);
        if( endptr == NULL ) return(false);
        if( *endptr == '\0' ) {
            value = lvalue;
            SetCurrentLineProcessed();
            return(true);
        }
    }

    return(false);
}

//------------------------------------------------------------------------------

bool CPrmFile::GetLogicalByKey(const CSmallString& name, bool& value)
{
    CSmallString svalue;

    // find key --------------------------------------
    if( FindKey(name,svalue) == true ) {
        if( svalue == NULL ) return(false);
        if( svalue == "on" || svalue == "ON" ||
            svalue == "true" || svalue == "TRUE" ||
            svalue == "T" || svalue == "1" ) {
            value = true;
            SetCurrentLineProcessed();
            return(true);
        }
        if( svalue == "off" || svalue == "OFF" ||
            svalue == "false" || svalue == "FALSE" ||
            svalue == "F" || svalue == "0" ) {
            value = false;
            SetCurrentLineProcessed();
            return(true);
        }
    }

    return(false);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================

bool CPrmFile::FindKey(const CSmallString& name, CSmallString& value)
{
    if( CurrentSection == NULL ) return(false);

    // find line with key
    CurrentLine = CurrentSection->FirstLine;

    CSmallString uname,lname,lvalue;

    uname = name;
    uname.ToUpperCase();

    while( CurrentLine != NULL ) {
        if( CurrentLine->Split(lname,lvalue) == false ) {
            // if no key value, try only key;
            lname = CurrentLine->GetText();
            lvalue = "";
        }
        lname.ToUpperCase();
        if( lname == uname ) {
            value = lvalue;
            return(true);
        }
        CurrentLine = CurrentLine->Next;
    }

    return(false);
}

//==============================================================================
//------------------------------------------------------------------------------
//==============================================================================


