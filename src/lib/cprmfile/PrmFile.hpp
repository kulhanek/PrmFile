#ifndef PrmFileH
#define PrmFileH
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

#include <SmallString.hpp>
#include <string>

//------------------------------------------------------------------------------

class CPrmGroup;
class CPrmSection;
class CPrmLine;

//------------------------------------------------------------------------------

class CPrmFile {
public:
// constructor/destructor -----------------------------------------------------
    //! constructor
    CPrmFile(void);
    //! destructor
    virtual ~CPrmFile(void);

// general prmfile methods ----------------------------------------------------
    //! discard previous data
    void Clear(void);

    //! read file - previous contents is deleted
    bool Read(const CSmallString& name);

    //! parse file
    bool Parse(FILE* p_file);

    //! print file - or print unprocessed items
    bool Dump(FILE* fout=NULL,bool unprocessed=false);

    //! print file - or print unprocessed items
    bool Dump(std::ostream& out,bool unprocessed=false);

    //! print group 'gname' - or print unprocessed items in it
    bool DumpGroup(const CSmallString& gname,FILE* fout=NULL,bool unprocessed=false);

    //! print group 'gname' - or print unprocessed items in it
    bool DumpGroup(const CSmallString& gname,std::ostream& out,bool unprocessed=false);

    //! count unprocessed lines in file
    int CountULines(void);

    //! count unprocessed lines in group 'gname'
    int CountULines(const CSmallString& gname);

// group methods --------------------------------------------------------------
    //! rewind to the first group
    bool FirstGroup(void);

    //! open group with name 'gname'
    bool OpenGroup(const CSmallString& gname);

    //! move to next group
    bool NextGroup(void);

    //! return number of section in the group
    int CountGroup(void);

    //! return name of current group
    const CSmallString& GetGroupName(void);

// section methods ------------------------------------------------------------
    //! rewind to the first section
    bool FirstSection(void);

    //! open section with name 'sname'
    bool OpenSection(const CSmallString& sname);

    //! open section with name 'sname' and return number of items
    int OpenSectionAndCount(const CSmallString& sname);

    //! move to next section
    bool NextSection(void);

    //! return number of items in the section
    int CountSection(void);

    //! return name of current section
    const CSmallString& GetSectionName(void);

    //! set whole section as processed
    void SetSecAsProcessed(void);

// line methods ---------------------------------------------------------------
    //! rewind to the first line of section
    bool FirstLine(void);

    //! move to next line
    bool NextLine(void);

    //! return key from current line
    bool GetLineKey(CSmallString& key);

    //! return value from current line
    bool GetLineValue(CSmallString& value);

    //! return the whole line
    bool GetLine(CSmallString& value);

    //! return the whole line
    bool GetLine(std::string& value);

    //! set current line as processed
    void SetCurrentLineProcessed(void);

// get values by key ----------------------------------------------------------
    bool GetIntegerByKey(const CSmallString& name, int& value);
    bool GetStringByKey(const CSmallString& name, CSmallString& value);
    bool GetStringByKey(const std::string& name, std::string& value);
    bool GetFloatByKey(const CSmallString& name, float& value);
    bool GetDoubleByKey(const CSmallString& name, double& value);
    bool GetLogicalByKey(const CSmallString& name, bool& value);

// section of private data ----------------------------------------------------
private:
    CSmallString    Name;               // name of file
    CPrmGroup*      _FirstGroup;        // first group in file
    CPrmGroup*      CurrentGroup;       // current group
    CPrmSection*    CurrentSection;     // current section
    CPrmLine*       CurrentLine;        // current line
    int             FieldPosition;      // used by prmfile_get_field

    //! parse line
    bool ParseLine(int lineid,const CSmallString& sline);

    //! add new group to a list
    bool AddGroup(int lineid,const CSmallString& gname);

    //! add new section to a list
    bool AddSection(int lineid,const CSmallString& sname);

    //! add line to a list
    bool AddLine(int lineid,const CSmallString& line);

    //! print error occured during file parsing
    void AddError(int lineid,const CSmallString& err);

    //! print group - or print unprocessed items in it
    bool DumpGroup(CPrmGroup* p_grp,std::ostream& out,const bool unprocessed);

    //! find key
    bool FindKey(const CSmallString& name, CSmallString& value);

    // static data
    static const char* white_list;
    static const char* comm_list;
    static const char* quot_list;

    friend class CPrmLine;
};

//------------------------------------------------------------------------------
#endif
